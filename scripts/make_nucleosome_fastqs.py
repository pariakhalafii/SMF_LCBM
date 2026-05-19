"""Generate a 550 bp synthetic SMF dataset with three positioned nucleosomes.

This is the "textbook" SMF picture you'd expect to see on a designed construct:

    [linker][----- nuc 1 -----][linker][----- nuc 2 -----][linker][----- nuc 3 -----][linker]
       0--26 |     27..147     |148-193|    194..341     |342-360|    361..508     |509-549

The biology being simulated:

* The methyltransferase **M.CviPI** can only reach **accessible** GpC sites, so it
  methylates GpCs in the **linker** DNA.
* GpC sites wrapped around histone octamers (the nucleosome regions above) are
  **inaccessible** and remain unmethylated.
* After bisulfite/EM-seq conversion, **unmethylated C -> T** while methylated C stays C.
* Net result on a read covering this construct:

    - linker GpCs:  preserved as **C** (~95% of the time, with realistic noise)
    - nucleosome GpCs: read as **T** (~95% of the time)
    - all other Cs (non-GpC contexts): not methylated by M.CviPI so they convert to T

So when you align the reads and look in IGV, you should see C's preserved at GpC
positions in the linker stretches, and T's everywhere else.

Outputs (in --out-dir, default data/nucleosomes/):

    reference.fa
    nuc_R1.fastq.gz
    nuc_R2.fastq.gz
    truth.tsv         <- per-read insert coordinates and strand, for cross-checking

Usage:

    python scripts/make_nucleosome_fastqs.py --n-reads 500 --seed 0
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
# The construct: 3 positioned nucleosomes inside a 550 bp molecule.
# Coordinates below are 1-based, inclusive on both ends, matching how the user
# described the experiment. We convert to 0-based half-open internally.
# ---------------------------------------------------------------------------

CONSTRUCT_LENGTH = 550

# Bismark needs a few bp of flanking reference around each read end to determine
# the trinucleotide context of every cytosine. If reads cover the entire reference
# exactly, Bismark discards them with "genomic sequence could not be extracted".
# So we wrap the 550 bp construct in flanking sequence to make a longer reference.
FLANK_BP = 75
REF_LENGTH = CONSTRUCT_LENGTH + 2 * FLANK_BP   # 700 bp total
CONSTRUCT_START_0 = FLANK_BP                    # 0-based start of the construct
CONSTRUCT_END_0 = CONSTRUCT_START_0 + CONSTRUCT_LENGTH  # exclusive

# Nucleosome coordinates as the user described (1-based, RELATIVE TO THE 550 BP
# CONSTRUCT). We shift by FLANK_BP to convert into reference coordinates.
NUCLEOSOMES_1BASED_CONSTRUCT = [
    ("nuc1",  27, 147),
    ("nuc2", 194, 341),
    ("nuc3", 361, 508),
]
NUCLEOSOMES_1BASED = [
    (name, s + FLANK_BP, e + FLANK_BP) for name, s, e in NUCLEOSOMES_1BASED_CONSTRUCT
]


def nucleosome_intervals_0based() -> list[tuple[str, int, int]]:
    """Return [(name, start, end_exclusive)] for the three nucleosomes in REFERENCE
    coordinates (0-based). The construct's nucleosome at 1-based 27-147 lives at
    reference 0-based positions (FLANK_BP+26, FLANK_BP+147)."""
    return [(name, s - 1, e) for name, s, e in NUCLEOSOMES_1BASED]


def is_in_nucleosome(pos: int) -> bool:
    """Is the 0-based REFERENCE position inside any of the three nucleosomes?"""
    for _, s, e in nucleosome_intervals_0based():
        if s <= pos < e:
            return True
    return False


# ---------------------------------------------------------------------------
# Reference construction -- random sequence sprinkled with GpCs every ~9 bp so
# the heatmap has nicely spaced columns. We embed real GpC dinucleotides without
# making the rest of the sequence GC-rich.
# ---------------------------------------------------------------------------

def make_reference(rng: np.random.Generator, length: int = REF_LENGTH,
                   gpc_spacing_bp: int = 9) -> str:
    """Random reference with all four bases (A/T-leaning), so that after bisulfite
    conversion the sequence still has enough complexity for bowtie2 to align.

    A pure-AT background (which we tried first) becomes essentially mononucleotide
    after bisulfite conversion and bowtie2 fails to align any reads.
    """
    # Slightly AT-leaning (matches typical mammalian genomes) but full complexity.
    seq = list(rng.choice(list("ACGT"), size=length, p=[0.30, 0.30, 0.20, 0.20]))

    # Force deliberate GpC sites at regular intervals so the SMF readout has
    # uniform coverage. Two ambiguous contexts must be avoided at every embedded site:
    #
    #   +strand GCG: would make the + strand C ambiguous (excluded by Bismark).
    #     -> force seq[i+2] != G
    #   -strand GCG (i.e., +strand CGC): would make the - strand call at this duplex
    #     site ambiguous, so only + strand reads contribute calls -> half-empty column
    #     in the per-read matrix -> visible white blocks in the heatmap.
    #     -> force seq[i-1] != C
    for i in range(3, length - 3, gpc_spacing_bp):
        # First, fix upstream and downstream neighbours BEFORE writing the GpC,
        # so we don't accidentally clobber a GpC we placed in a previous iteration.
        if seq[i - 1] == "C":
            seq[i - 1] = rng.choice(list("ATG"))   # non-C upstream -> avoid +strand CGC
        seq[i] = "G"
        seq[i + 1] = "C"
        if seq[i + 2] == "G":
            seq[i + 2] = rng.choice(list("ATC"))   # non-G downstream -> avoid +strand GCG

    # Iteratively scrub two kinds of ambiguous contexts in the RANDOM background:
    #   * +strand GCG  -> fix by mutating the trailing G to A/T
    #   * +strand CGC  -> equivalent to -strand GCG; fix by mutating the leading C to A/T
    # Use only A or T as replacements so we don't accidentally create a new G or C
    # adjacent to existing patterns. Iterate until the sequence is clean.
    for _ in range(20):
        changed = False
        for i in range(1, length - 1):
            if seq[i - 1] == "G" and seq[i] == "C" and seq[i + 1] == "G":
                seq[i + 1] = rng.choice(list("AT"))
                changed = True
            elif seq[i - 1] == "C" and seq[i] == "G" and seq[i + 1] == "C":
                # Don't touch the C if it's the C of an embedded GpC (then seq[i-2]==G);
                # in that case break the upstream side instead.
                if i >= 2 and seq[i - 2] == "G":
                    seq[i + 1] = rng.choice(list("AT"))
                else:
                    seq[i - 1] = rng.choice(list("AT"))
                changed = True
        if not changed:
            break

    return "".join(seq).upper()


# ---------------------------------------------------------------------------
# Methylation per molecule.
# Linker GpC: methylated with prob p_meth_linker     (default 0.95)
# Nucleosome GpC: methylated with prob p_meth_nuc    (default 0.05)
# Non-GpC Cs: never methylated by M.CviPI -- they always convert to T.
# ---------------------------------------------------------------------------

def gpc_positions(seq: str) -> list[int]:
    """Indices of every GpC site (return the + strand C index)."""
    return [i + 1 for i in range(len(seq) - 1) if seq[i] == "G" and seq[i + 1] == "C"]


def cpg_positions(seq: str) -> list[int]:
    """Indices of every CpG site (return the + strand C index)."""
    return [i for i in range(len(seq) - 1) if seq[i] == "C" and seq[i + 1] == "G"]


def decide_methylated(rng: np.random.Generator, seq: str,
                      p_meth_linker: float, p_meth_nuc: float) -> set[int]:
    """Return the set of + strand C positions methylated on this molecule.

    Dual-enzyme dSMF: M.CviPI methylates accessible GpC sites, M.SssI methylates
    accessible CpG sites. Both enzymes follow the same accessibility logic --
    linker DNA gets methylated, nucleosome-protected DNA does not.
    """
    out: set[int] = set()
    # M.CviPI -- accessibility readout at GpC sites
    for pos in gpc_positions(seq):
        p = p_meth_nuc if is_in_nucleosome(pos) else p_meth_linker
        if rng.random() < p:
            out.add(pos)
    # M.SssI -- accessibility readout at CpG sites
    for pos in cpg_positions(seq):
        p = p_meth_nuc if is_in_nucleosome(pos) else p_meth_linker
        if rng.random() < p:
            out.add(pos)
    return out


def bisulfite_convert(seq: str, methylated: set[int], strand: str) -> str:
    """Apply C->T to every cytosine that is NOT in `methylated`, on the given strand.

    For "+" we work in reference frame; for "-" we reverse-complement first and
    translate methylated positions accordingly.
    """
    if strand == "+":
        return "".join("T" if (b == "C" and i not in methylated) else b
                       for i, b in enumerate(seq))

    rc = revcomp(seq)
    n = len(seq)
    # Translate each + strand methylated C to its - strand counterpart on the rc
    # sequence. The mapping depends on context:
    #   * GpC (G at i-1, C at i):   - strand C of same duplex is at + strand pos i-1,
    #                                which lands at rc position n-1-(i-1) = n-i.
    #   * CpG (C at i, G at i+1):   - strand C of same duplex is at + strand pos i+1,
    #                                which lands at rc position n-1-(i+1) = n-i-2.
    rc_meth: set[int] = set()
    for i in methylated:
        if i > 0 and seq[i - 1] == "G":
            rc_meth.add(n - i)             # GpC mapping
        elif i < n - 1 and seq[i + 1] == "G":
            rc_meth.add(n - i - 2)         # CpG mapping
        # else: shouldn't happen since methylated set only contains GpC/CpG Cs
    return "".join("T" if (b == "C" and i not in rc_meth) else b
                   for i, b in enumerate(rc))


# ---------------------------------------------------------------------------
# FASTQ writing helpers
# ---------------------------------------------------------------------------

def phred_string(seq_len: int, rng: np.random.Generator, mean_q: int = 35) -> str:
    qs = rng.integers(mean_q - 3, mean_q + 3, size=seq_len)
    qs = np.clip(qs, 2, 41)
    return "".join(chr(int(q) + 33) for q in qs)


def write_fastq_record(fh, header: str, seq: str, qual: str) -> None:
    fh.write(f"@{header}\n{seq}\n+\n{qual}\n".encode())


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(argv: list[str]) -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out-dir", default="data/nucleosomes", help="output directory")
    p.add_argument("--n-reads", type=int, default=500)
    p.add_argument("--read-len", type=int, default=300)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--p-meth-linker", type=float, default=0.95,
                   help="P(methylation) at linker GpCs (default 0.95)")
    p.add_argument("--p-meth-nuc", type=float, default=0.05,
                   help="P(methylation) at nucleosome GpCs (default 0.05)")
    p.add_argument("--gpc-spacing-bp", type=int, default=9,
                   help="approximate spacing between embedded GpC sites (default 9)")
    args = p.parse_args(argv)

    if args.read_len * 2 < CONSTRUCT_LENGTH + 10:
        print(f"ERROR: 2 * read_len ({2*args.read_len}) must exceed construct length "
              f"({CONSTRUCT_LENGTH}) so R1/R2 overlap.", file=sys.stderr)
        return 2

    rng = np.random.default_rng(args.seed)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # 1. Reference: the 550 bp construct wrapped in flanks. The flanks give Bismark
    #    enough surrounding sequence to extract trinucleotide contexts at read edges.
    reference = make_reference(rng, length=REF_LENGTH, gpc_spacing_bp=args.gpc_spacing_bp)
    fa = out_dir / "reference.fa"
    with fa.open("w") as fh:
        fh.write(">nucleosome_construct\n")
        for i in range(0, len(reference), 70):
            fh.write(reference[i : i + 70] + "\n")

    n_gpc = len(gpc_positions(reference))
    n_gpc_nuc = sum(1 for p in gpc_positions(reference) if is_in_nucleosome(p))
    print(f"[nuc] reference: {fa}  ({len(reference)} bp"
          f" = {FLANK_BP} bp flank + {CONSTRUCT_LENGTH} bp construct + {FLANK_BP} bp flank)")
    print(f"[nuc] embedded GpC sites: {n_gpc} total, {n_gpc_nuc} inside nucleosomes, "
          f"{n_gpc - n_gpc_nuc} elsewhere")
    print(f"[nuc] nucleosome positions in REFERENCE coordinates (1-based, inclusive):")
    for name, s, e in NUCLEOSOMES_1BASED:
        print(f"[nuc]   {name}: {s}-{e}  ({e - s + 1} bp)")
    print(f"[nuc]   (the user's construct-relative positions 27-147, 194-341, 361-508")
    print(f"[nuc]    are these +{FLANK_BP} after wrapping in flanks)")

    # 2. Reads -- each molecule is the 550 bp construct (positions FLANK_BP..FLANK_BP+550
    #    of the reference), with bisulfite/EM-seq conversion applied per molecule.
    strands = rng.choice(["+", "-"], size=args.n_reads)
    overlap = 2 * args.read_len - CONSTRUCT_LENGTH

    r1_path = out_dir / "nuc_R1.fastq.gz"
    r2_path = out_dir / "nuc_R2.fastq.gz"
    truth_path = out_dir / "truth.tsv"

    with gzip.open(r1_path, "wb") as r1_fh, \
         gzip.open(r2_path, "wb") as r2_fh, \
         truth_path.open("w") as truth_fh:

        truth_fh.write("read_id\tstrand\tinsert_start\tinsert_end\tn_meth_gpc\tn_unmeth_gpc\n")

        for i in range(args.n_reads):
            strand = str(strands[i])
            methylated = decide_methylated(rng, reference,
                                           p_meth_linker=args.p_meth_linker,
                                           p_meth_nuc=args.p_meth_nuc)
            converted = bisulfite_convert(reference, methylated, strand)

            # Extract just the construct portion as the molecule.
            if strand == "+":
                molecule = converted[CONSTRUCT_START_0:CONSTRUCT_END_0]
            else:
                # The converted "-" strand is the rc of the FULL reference. The
                # construct on the - strand corresponds to the same coordinates but
                # measured from the - strand 5' end.
                rc_construct_start = REF_LENGTH - CONSTRUCT_END_0
                rc_construct_end = REF_LENGTH - CONSTRUCT_START_0
                molecule = converted[rc_construct_start:rc_construct_end]

            r1 = molecule[: args.read_len]
            r2 = revcomp(molecule[-args.read_len:])

            read_id = f"NUC:1:1:{i+1}:{i+1}"
            write_fastq_record(r1_fh, f"{read_id} 1:N:0:1", r1, phred_string(len(r1), rng))
            write_fastq_record(r2_fh, f"{read_id} 2:N:0:1", r2, phred_string(len(r2), rng))

            n_meth = sum(1 for p in methylated if not is_in_nucleosome(p))  # outside nucleosomes
            n_unmeth = sum(1 for p in gpc_positions(reference)
                           if is_in_nucleosome(p) and p not in methylated)
            truth_fh.write(
                f"{read_id}\t{strand}\t{CONSTRUCT_START_0}\t{CONSTRUCT_END_0}\t{n_meth}\t{n_unmeth}\n"
            )

    print(f"[nuc] wrote {args.n_reads} read pairs:")
    print(f"        R1: {r1_path}  ({args.read_len} bp each)")
    print(f"        R2: {r2_path}  ({args.read_len} bp each)")
    print(f"     truth: {truth_path}")
    print(f"[nuc] R1/R2 overlap: ~{overlap} bp on a {CONSTRUCT_LENGTH} bp insert "
          f"(FLASH min-overlap=30 will catch this)")
    print()
    print("[nuc] Next steps:")
    print(f"    1. Rebuild the Bismark index (the reference changed):")
    print(f"         rm -rf {fa.parent}/Bisulfite_Genome")
    print(f"         bismark_genome_preparation --bowtie2 {fa.parent}/")
    print(f"    2. Drive the pipeline with config/config.nucleosomes.yaml")
    print(f"    3. Use --feature-center {NUCLEOSOMES_1BASED[0][1] + (NUCLEOSOMES_1BASED[0][2]-NUCLEOSOMES_1BASED[0][1])//2} "
          f"(nuc1), {NUCLEOSOMES_1BASED[1][1] + (NUCLEOSOMES_1BASED[1][2]-NUCLEOSOMES_1BASED[1][1])//2} "
          f"(nuc2), {NUCLEOSOMES_1BASED[2][1] + (NUCLEOSOMES_1BASED[2][2]-NUCLEOSOMES_1BASED[2][1])//2} "
          f"(nuc3) for the analysis.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
