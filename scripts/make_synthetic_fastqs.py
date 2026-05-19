"""Generate a tiny synthetic SMF dataset (reference + paired-end FASTQs).

The point of this script is to give you something to feed the actual upstream
pipeline (FLASH -> Bismark -> samtools -> custom caller) WITHOUT using your real
data. It produces a few hundred reads in a few seconds, all formatted exactly
like real Illumina output, so the shell scripts run unchanged.

What it builds:

  1. A fake reference (~620 bp) with two pretend TF motifs and lots of GpC sites.
  2. N synthetic molecules, each labelled with one of:
        accessible  -- every GpC methylated
        tf_bound    -- GpCs methylated *except* a small patch over the motif
        nucleosome  -- GpCs in a wide ~150 bp patch all unmethylated
  3. Bisulfite/EM-seq conversion: unmethylated C -> T on the chosen strand.
  4. Each converted molecule is split into:
        R1 = first 300 bp
        R2 = reverse complement of the last 300 bp
     With a 550 bp insert, R1 and R2 overlap by ~50 bp (matches real 2x300
     Illumina on 3N601-style amplicons; FLASH's min-overlap=30 will catch them).
  5. Phred-35-ish quality scores (real-looking but cleaner than real data).

Outputs (in --out-dir, default data/synthetic/):
    reference.fa
    synth_R1.fastq.gz
    synth_R2.fastq.gz
    truth.tsv           <- one row per read with the ground-truth state

Then point config/config.synthetic.yaml at these and run the pipeline.

Usage:
    python scripts/make_synthetic_fastqs.py --n-reads 500 --seed 0
"""
from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from smf.context import revcomp  # noqa: E402


# ---------------------------------------------------------------------------
# Reference construction
# ---------------------------------------------------------------------------

# Two arbitrary "TF motifs" we'll embed in the reference so the classifier has
# something specific to find. Real protocols would use actual motifs (CTCF, etc.).
MOTIF_A = "TGACTCA"  # at position ~280
MOTIF_B = "GATAAG"   # at position ~430


def make_reference(rng: np.random.Generator, length: int = 620) -> tuple[str, dict]:
    """Build a fake reference sequence with embedded motifs and lots of GpCs."""
    # Start with random A/T-heavy sequence, then sprinkle GC-rich windows so we
    # get plenty of informative GpC sites.
    base_pool = list("ATATATATATGCGCGCGC")  # ~22% G/C overall
    seq = list(rng.choice(base_pool, size=length))

    # Sprinkle in extra GpC sites every ~10 bp so the heatmap has dense coverage.
    for i in range(5, length - 2, 10):
        seq[i] = "G"
        seq[i + 1] = "C"
        # Add a non-G base after to keep it as GCH (informative GpC, not GCG).
        seq[i + 2] = rng.choice(list("ATC"))

    # Embed the motifs.
    pos_a = 280
    pos_b = 430
    for j, b in enumerate(MOTIF_A):
        seq[pos_a + j] = b
    for j, b in enumerate(MOTIF_B):
        seq[pos_b + j] = b

    full = "".join(seq).upper()
    motif_centers = {
        "motif_A": pos_a + len(MOTIF_A) // 2,
        "motif_B": pos_b + len(MOTIF_B) // 2,
    }
    return full, motif_centers


# ---------------------------------------------------------------------------
# Methylation and conversion
# ---------------------------------------------------------------------------

def _gpc_positions(seq: str) -> list[int]:
    """Indices of every GpC site (G followed by C, treating the C as the call site).

    For SMF the *C* of a GpC carries the methylation, so we return that C index.
    """
    out = []
    for i in range(len(seq) - 1):
        if seq[i] == "G" and seq[i + 1] == "C":
            out.append(i + 1)
    return out


def _is_protected(pos: int, label: str, centers: dict, cfg) -> bool:
    """Decide if a GpC at `pos` is in the protected footprint for this state."""
    if label == "accessible":
        return False
    if label == "tf_bound":
        # Small patch over either motif.
        for c in centers.values():
            if abs(pos - c) <= cfg.tf_half_width:
                return True
        return False
    if label == "nucleosome":
        # Wide patch over motif A (one nucleosome's worth).
        return abs(pos - centers["motif_A"]) <= cfg.nuc_half_width
    return False


def _bisulfite_convert(seq: str, methylated: set[int], strand: str) -> str:
    """Apply C->T conversion to unmethylated cytosines on the given strand."""
    if strand == "+":
        out = []
        for i, b in enumerate(seq):
            if b == "C" and i not in methylated:
                out.append("T")
            else:
                out.append(b)
        return "".join(out)
    # On the bottom strand we flip first, convert, leave the read in that frame.
    rc = revcomp(seq)
    n = len(seq)
    # Translate plus-strand methylated set into bottom-strand coordinates.
    rc_meth = {n - 1 - i for i in methylated}
    out = []
    for i, b in enumerate(rc):
        if b == "C" and i not in rc_meth:
            out.append("T")
        else:
            out.append(b)
    return "".join(out)


# ---------------------------------------------------------------------------
# FASTQ writing
# ---------------------------------------------------------------------------

def _phred(seq_len: int, rng: np.random.Generator, mean_q: int = 35) -> str:
    """A realistic-but-clean Phred quality string."""
    qs = rng.integers(mean_q - 3, mean_q + 3, size=seq_len)
    qs = np.clip(qs, 2, 41)
    return "".join(chr(int(q) + 33) for q in qs)


def _write_fastq_record(fh, read_id: str, seq: str, qual: str) -> None:
    fh.write(f"@{read_id}\n{seq}\n+\n{qual}\n".encode())


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

class _Cfg:
    def __init__(self, args):
        self.n_reads = args.n_reads
        self.read_len = args.read_len
        self.insert_len = args.insert_len
        self.tf_half_width = 12
        self.nuc_half_width = 75
        self.frac_accessible = 0.4
        self.frac_tf_bound = 0.25
        self.frac_nucleosome = 0.3
        self.p_meth_accessible = 0.9
        self.p_meth_protected = 0.05


def main(argv: list[str]) -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out-dir", default="data/synthetic", help="output directory")
    p.add_argument("--n-reads", type=int, default=500, help="number of read pairs to generate")
    p.add_argument("--read-len", type=int, default=300, help="length of R1 and R2 (bp)")
    p.add_argument("--insert-len", type=int, default=550,
                   help="length of the merged molecule (must be < 2*read_len so R1/R2 overlap)")
    p.add_argument("--seed", type=int, default=0)
    args = p.parse_args(argv)

    if args.insert_len >= 2 * args.read_len:
        print("ERROR: insert_len must be < 2 * read_len so R1/R2 overlap", file=sys.stderr)
        return 2

    cfg = _Cfg(args)
    rng = np.random.default_rng(args.seed)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # 1. Reference
    ref_len = args.insert_len + 70   # a little extra so molecules can sit at varying positions
    reference, centers = make_reference(rng, length=ref_len)
    fa = out_dir / "reference.fa"
    with fa.open("w") as fh:
        fh.write(">synth_construct\n")
        for i in range(0, len(reference), 70):
            fh.write(reference[i : i + 70] + "\n")
    print(f"[synth] reference: {fa}  ({len(reference)} bp)")
    print(f"[synth] motif centres: {centers}")

    # 2. Molecules
    fracs = np.array([cfg.frac_accessible, cfg.frac_tf_bound, cfg.frac_nucleosome])
    fracs = np.append(fracs, max(0.0, 1.0 - fracs.sum()))
    fracs = fracs / fracs.sum()
    labels = rng.choice(
        ["accessible", "tf_bound", "nucleosome", "unclassified"],
        size=cfg.n_reads, p=fracs,
    )
    strands = rng.choice(["+", "-"], size=cfg.n_reads)
    gpc_sites = _gpc_positions(reference)

    r1_path = out_dir / "synth_R1.fastq.gz"
    r2_path = out_dir / "synth_R2.fastq.gz"
    truth_path = out_dir / "truth.tsv"

    with gzip.open(r1_path, "wb") as r1_fh, \
         gzip.open(r2_path, "wb") as r2_fh, \
         truth_path.open("w") as truth_fh:

        truth_fh.write("read_id\tstate\tstrand\tinsert_start\tinsert_end\n")

        for i in range(cfg.n_reads):
            label = str(labels[i])
            strand = str(strands[i])

            # Decide methylation per GpC site for this molecule.
            methylated: set[int] = set()
            for pos in gpc_sites:
                p_meth = (cfg.p_meth_protected
                          if _is_protected(pos, label, centers, cfg)
                          else cfg.p_meth_accessible)
                if rng.random() < p_meth:
                    methylated.add(pos)

            # Pick where in the reference this molecule sits.
            max_start = ref_len - cfg.insert_len
            insert_start = int(rng.integers(0, max_start + 1))
            insert_end = insert_start + cfg.insert_len
            insert_seq = reference[insert_start:insert_end]

            # Translate methylated positions into insert-local coordinates.
            local_meth = {p - insert_start for p in methylated
                          if insert_start <= p < insert_end}

            converted = _bisulfite_convert(insert_seq, local_meth, strand)

            # Build R1 and R2 from the converted insert.
            r1 = converted[: cfg.read_len]
            r2 = revcomp(converted[-cfg.read_len:])

            read_id = f"SYNTH:1:1:{i+1}:{i+1}"
            _write_fastq_record(r1_fh, f"{read_id} 1:N:0:1", r1, _phred(len(r1), rng))
            _write_fastq_record(r2_fh, f"{read_id} 2:N:0:1", r2, _phred(len(r2), rng))

            truth_fh.write(f"{read_id}\t{label}\t{strand}\t{insert_start}\t{insert_end}\n")

    print(f"[synth] wrote {cfg.n_reads} read pairs:")
    print(f"        R1: {r1_path}")
    print(f"        R2: {r2_path}")
    print(f"     truth: {truth_path}")
    print()
    print("[synth] Next steps:")
    print(f"    1. Build a Bismark genome index for {fa.parent}/")
    print(f"       (the prep tool expects a *folder* containing the FASTA):")
    print(f"         bismark_genome_preparation --bowtie2 {fa.parent}/")
    print(f"    2. Use config/config.synthetic.yaml to drive the pipeline.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
