"""Run the SMF downstream analysis on a per-read calls TSV.

This script picks up where ``scripts/extract_per_read_methylation.py`` leaves off:
takes the TSV of per-read methylation calls, builds a per-read matrix, classifies
single-molecule footprint states, and writes plots + summaries to disk.

If you point ``--truth`` at a ground-truth TSV (e.g. the one produced by
``make_synthetic_fastqs.py``), it also prints a confusion matrix so you can see how
well the classifier recovers the simulated states.

Usage:

    python scripts/analyze_per_read_calls.py \\
        --tsv     results/synthetic/methylation/per_read_calls.tsv \\
        --feature-center 283 \\
        --out-dir results/synthetic/plots \\
        --truth   data/synthetic/truth.tsv
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from smf import classify, io as smf_io, per_read, viz  # noqa: E402


def _strip_pair_suffix(read_id: str) -> str:
    """Bismark sometimes appends `_1:N:0:1` (the FASTQ comment) to the read name.
    Strip that off so we can match against truth files keyed by the bare read ID."""
    return re.sub(r"_[12]:N:\d+:\d+$", "", read_id)


def main(argv: list[str]) -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--tsv", required=True, help="per-read calls TSV from extract_per_read_methylation.py")
    p.add_argument("--feature-center", type=int, required=True,
                   help="reference position (1-based) of the feature you want to classify around")
    p.add_argument("--contexts", default="GCH,HCG",
                   help="comma-separated contexts to use. Default GCH,HCG = both readouts "
                        "from dual-enzyme dSMF (M.CviPI on GpC + M.SssI on CpG). "
                        "Pass just 'GCH' for GpC-only accessibility, just 'HCG' for CpG-only.")
    p.add_argument("--sort-by", choices=["feature", "global", "none"], default="global",
                   help="how to sort the rows of the heatmap. "
                        "'feature' sorts by accessibility in a window around --feature-center "
                        "(differs across feature centers). "
                        "'global' (default) sorts by mean accessibility across the whole window "
                        "(SAME row order regardless of feature center, useful for comparing plots). "
                        "'none' keeps the input order.")
    p.add_argument("--out-dir", required=True, help="where to write plots + CSV summaries")
    p.add_argument("--truth", default=None,
                   help="optional ground-truth TSV with read_id and state columns")
    p.add_argument("--free-dna-threshold", type=float, default=0.65,
                   help="mean GpC methylation in the nucleosome zone above which a read is "
                        "classified as free DNA (default 0.65). Set to 1.0 to disable filtering.")
    p.add_argument("--nuc-half-width", type=int, default=73,
                   help="half-width in bp of the nucleosome zone used for free-DNA detection "
                        "(default 73 bp, i.e. ±73 bp around --feature-center).")
    args = p.parse_args(argv)

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    contexts = tuple(c.strip() for c in args.contexts.split(",") if c.strip())

    # -------------------------------------------------- 1. Load
    print(f"[analyze] loading {args.tsv}")
    prm = smf_io.read_per_read_tsv(args.tsv, contexts=contexts)
    print(f"  {prm.n_reads} reads x {prm.n_sites} informative {','.join(contexts)} sites")

    # Build a per-site context lookup for the coverage plot colouring.
    _raw = pd.read_csv(args.tsv, sep="\t", usecols=["pos", "context"],
                       dtype={"context": str})
    _site_ctx = _raw[_raw["context"].isin(contexts)].groupby("pos")["context"].first()
    site_contexts = [_site_ctx.get(int(s), "GCH") for s in prm.sites]

    # -------------------------------------------------- 2. (no coverage filter)
    prm_f = prm

    # -------------------------------------------------- 3. Sort + classify
    if args.sort_by == "feature":
        sorted_prm = per_read.sort_by_pattern(prm_f, feature_center=args.feature_center, window=80)
        sort_desc = f"sorted by accessibility around {args.feature_center}"
    elif args.sort_by == "global":
        # Sort by mean accessibility across the whole window. This produces the SAME
        # row order regardless of which feature centre is being analyzed, so plots
        # for different features over the same data are directly comparable.
        sorted_prm = per_read.sort_by_pattern(
            prm_f,
            feature_center=int((prm_f.region_start + prm_f.region_end) / 2),
            window=int((prm_f.region_end - prm_f.region_start) / 2 + 1),
        )
        sort_desc = "sorted globally (same row order across feature centres)"
    else:
        sorted_prm = prm_f
        sort_desc = "input order"

    # -------------------------------------------------- 3b. Free-DNA detection
    # Compute mask in prm_f order (for filtering) and sorted_prm order (for flagged heatmap).
    free_dna_mask = prm_f.compute_free_dna_mask(
        feature_center=args.feature_center,
        nuc_half_width=args.nuc_half_width,
        threshold=args.free_dna_threshold,
    )
    free_dna_mask_sorted = sorted_prm.compute_free_dna_mask(
        feature_center=args.feature_center,
        nuc_half_width=args.nuc_half_width,
        threshold=args.free_dna_threshold,
    )
    n_free = int(free_dna_mask.sum())
    n_total = prm_f.n_reads
    print(f"[analyze] free-DNA filter: {n_free}/{n_total} reads flagged "
          f"(mean GpC in ±{args.nuc_half_width} bp > {args.free_dna_threshold})")
    prm_clean = prm_f.filter_by_mask(~free_dna_mask)

    if args.sort_by == "feature":
        sorted_prm_clean = per_read.sort_by_pattern(prm_clean, feature_center=args.feature_center, window=80)
    elif args.sort_by == "global":
        sorted_prm_clean = per_read.sort_by_pattern(
            prm_clean,
            feature_center=int((prm_clean.region_start + prm_clean.region_end) / 2),
            window=int((prm_clean.region_end - prm_clean.region_start) / 2 + 1),
        )
    else:
        sorted_prm_clean = prm_clean

    result = classify.classify_footprints(prm_clean, feature_center=args.feature_center)
    print(f"[analyze] heatmap rows {sort_desc}")
    print("[analyze] classifier counts (free-DNA excluded):", result.counts())

    # -------------------------------------------------- 4. Plots
    print("[analyze] writing plots to", out_dir)

    # Heatmap 1: free DNA filtered out entirely
    viz.single_molecule_heatmap(
        sorted_prm_clean, feature_center=args.feature_center,
        title=f"SMF heatmap (free DNA removed, n={prm_clean.n_reads}) -- centred at {args.feature_center}",
    ).savefig(out_dir / "single_molecule_heatmap.png", dpi=150)

    # Heatmap 2: all reads, free DNA highlighted green in their natural sorted position
    viz.single_molecule_heatmap(
        sorted_prm, feature_center=args.feature_center,
        title=f"SMF heatmap (free DNA highlighted, n={prm_f.n_reads}) -- centred at {args.feature_center}",
        free_dna_mask=free_dna_mask_sorted,
        sort_free_dna=False,
    ).savefig(out_dir / "single_molecule_heatmap_flagged.png", dpi=150)

    viz.average_accessibility_plot(
        prm_f, feature_center=args.feature_center, smooth_bp=20,
        title="Average GpC accessibility",
    ).savefig(out_dir / "average_accessibility.png", dpi=150)

    viz.coverage_per_position_plot(
        prm_f, feature_center=args.feature_center, low_coverage_threshold=20,
        site_contexts=site_contexts,
    ).savefig(out_dir / "coverage_per_position.png", dpi=150)

    viz.state_composition_plot(result.counts()).savefig(
        out_dir / "state_composition.png", dpi=150)

    viz.read_alignment_plot(
        prm_f, feature_center=args.feature_center,
        title=f"Read alignment (n={prm_f.n_reads}) -- centred at {args.feature_center}",
    ).savefig(out_dir / "read_alignment.png", dpi=150)

    # -------------------------------------------------- 5. Per-read summary
    summary = pd.DataFrame({
        "read_id": prm_clean.read_ids,
        "predicted_state": result.states,
        "motif_score": result.motif_score,
        "inner_flank_score": result.inner_flank_score,
        "outer_flank_score": result.outer_flank_score,
    })

    # -------------------------------------------------- 6. Compare to truth (optional)
    if args.truth:
        truth = pd.read_csv(args.truth, sep="\t")
        if "state" not in truth.columns or "read_id" not in truth.columns:
            print(f"[analyze] WARNING: --truth must have read_id and state columns; got {list(truth.columns)}")
        else:
            # Normalise read IDs on both sides (Bismark may have mangled the headers).
            truth["read_id_norm"] = truth["read_id"].astype(str).str.replace(" ", "_", regex=False)
            summary["read_id_norm"] = summary["read_id"].apply(_strip_pair_suffix)
            merged = summary.merge(truth[["read_id_norm", "state"]], on="read_id_norm", how="left")
            merged.rename(columns={"state": "true_state"}, inplace=True)

            confusion = pd.crosstab(
                merged["true_state"].fillna("missing"),
                merged["predicted_state"].fillna("missing"),
                dropna=False,
            )
            print("\n[analyze] confusion matrix (rows = truth, cols = prediction):")
            print(confusion)
            confusion.to_csv(out_dir / "confusion_matrix.csv")

            summary = merged.drop(columns=["read_id_norm"])

    summary.to_csv(out_dir / "per_read_summary.csv", index=False)
    print(f"[analyze] wrote per-read summary to {out_dir/'per_read_summary.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
