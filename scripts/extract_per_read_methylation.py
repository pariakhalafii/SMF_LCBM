"""Custom SMF per-read methylation caller.

Why we need this on top of Bismark's bulk caller:

* Bismark's ``bismark_methylation_extractor`` collapses everything to per-cytosine
  *averages* across all reads. For SMF we need the **per-read** binary call at every
  informative cytosine on the molecule -- that's what becomes a row in the heatmap.
* Bismark calls CpG, CHG, CHH but doesn't natively split out the **GpC** sites that
  carry the SMF accessibility signal. We resolve context here from the reference so
  GpC, CpG, GCG (ambiguous, dropped), and other-CHH/CHG are correctly partitioned.

Inputs
------
* sorted, indexed Bismark BAM (single-end on FLASH-merged reads)
* reference FASTA (same one used for ``bismark_genome_preparation``)

Output
------
A long-format TSV with one row per (read, informative cytosine):

    read_id  chrom  pos  strand  ref_base  context  call  methylated  read_length

Plus, if ``--molecule-pattern`` is given, a wide companion file with one row per read
containing the methylation pattern (a string of ``M`` / ``U`` / ``.`` characters at
every informative cytosine in genomic order).

Usage
-----

    python scripts/extract_per_read_methylation.py \\
        --bam   results/demo/align/*.minlen500.sorted.bam \\
        --fasta refs/genome/genome.fa \\
        --region chr1:10000-12000 \\
        --contexts GCH,HCG \\
        --out results/demo/methylation/per_read_calls.tsv \\
        --molecule-pattern results/demo/methylation/molecule_patterns.tsv
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from smf import io as smf_io  # noqa: E402
from smf.context import Context  # noqa: E402


def parse_region(region: str) -> tuple[str, int, int]:
    """Parse ``chr1:10000-12000`` -> ``("chr1", 9999, 12000)`` (0-based half-open)."""
    if ":" not in region:
        raise ValueError(f"region must be chrom:start-end, got {region!r}")
    chrom, span = region.split(":", 1)
    start_s, end_s = span.replace(",", "").split("-", 1)
    start = int(start_s) - 1   # input is 1-based inclusive (samtools convention)
    end = int(end_s)
    if end <= start:
        raise ValueError(f"region end <= start: {region}")
    return chrom, start, end


def main(argv: list[str]) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bam", required=True, help="sorted, indexed Bismark BAM")
    p.add_argument("--fasta", required=True, help="reference FASTA (with .fai index)")
    p.add_argument("--region", required=True, help="chrom:start-end (1-based inclusive)")
    p.add_argument(
        "--contexts",
        default="GCH",
        help="comma-separated contexts to keep (any of GCH,HCG,GCG,CCG,CHH,CHG). "
             "Default: GCH (the SMF accessibility readout). For endogenous CpG add HCG.",
    )
    p.add_argument("--min-mapq", type=int, default=10, help="minimum MAPQ (default 10)")
    p.add_argument("--out", required=True, help="output TSV with per-read calls")
    p.add_argument(
        "--molecule-pattern",
        default=None,
        help="optional output TSV with one row per read = pattern string of M/U/. chars",
    )
    args = p.parse_args(argv)

    chrom, start, end = parse_region(args.region)
    wanted = {Context(c.strip()) for c in args.contexts.split(",") if c.strip()}

    try:
        import pyfaidx  # type: ignore[import-untyped]
    except ImportError:
        print("ERROR: pyfaidx is required (pip install pyfaidx)", file=sys.stderr)
        return 2

    fasta = pyfaidx.Fasta(args.fasta, sequence_always_upper=True)

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    n_calls = 0
    by_read: dict[str, dict[int, str]] = {}  # read_id -> {pos: 'M'|'U'}
    read_meta: dict[str, dict[str, object]] = {}

    with out_path.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow([
            "read_id", "chrom", "pos", "strand", "ref_base",
            "context", "call", "methylated",
        ])

        for call in smf_io.iter_xm_calls(
            args.bam, region=(chrom, start, end), fasta=fasta, min_mapq=args.min_mapq,
        ):
            if call.context not in wanted:
                continue
            ref_base = "C" if call.strand == "+" else "G"
            writer.writerow([
                call.read_id, call.chrom, call.pos + 1, call.strand, ref_base,
                call.context.value,
                "methylated" if call.methylated else "unmethylated",
                int(call.methylated),
            ])
            n_calls += 1

            if args.molecule_pattern:
                by_read.setdefault(call.read_id, {})[call.pos] = "M" if call.methylated else "U"
                read_meta.setdefault(call.read_id, {"chrom": call.chrom, "strand": call.strand})

    print(f"[caller] wrote {n_calls} per-read calls to {out_path}")

    if args.molecule_pattern:
        all_positions = sorted({pos for d in by_read.values() for pos in d})
        pat_path = Path(args.molecule_pattern)
        pat_path.parent.mkdir(parents=True, exist_ok=True)
        with pat_path.open("w", newline="") as fh:
            writer = csv.writer(fh, delimiter="\t")
            header = ["read_id", "chrom", "strand", "n_calls"] + [str(p + 1) for p in all_positions]
            writer.writerow(header)
            for rid, calls in by_read.items():
                pattern = [calls.get(p, ".") for p in all_positions]
                meta = read_meta[rid]
                writer.writerow([rid, meta["chrom"], meta["strand"], sum(1 for c in pattern if c != "."), *pattern])
        print(f"[caller] wrote molecule patterns ({len(by_read)} reads) to {pat_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
