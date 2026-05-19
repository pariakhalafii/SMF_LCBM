"""End-to-end downstream SMF example, using only synthetic data.

Run from the project root:

    python scripts/run_example.py

Writes three plots and a summary CSV into ``examples/``.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

# Make the project importable when this script is run as `python scripts/run_example.py`.
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from scripts.make_synthetic_data import SyntheticConfig, make_synthetic  # noqa: E402
from smf import classify, per_read, stats, viz  # noqa: E402


def main() -> int:
    out_dir = ROOT / "examples"
    out_dir.mkdir(exist_ok=True)

    print("[example] Generating synthetic per-read GpC matrix")
    cfg = SyntheticConfig(seed=42)
    prm, true_labels, fragment_lengths = make_synthetic(cfg)
    print(f"  {prm.n_reads} reads x {prm.n_sites} sites, "
          f"window = [{prm.region_start}, {prm.region_end}], centre = {cfg.feature_center}")

    # --- Fragment-length QC (before filtering) ---
    fig_len = viz.fragment_length_plot(fragment_lengths, min_keep=180,
                                       title="Synthetic SMF -- fragment length distribution")
    fig_len.savefig(out_dir / "fragment_length.png", dpi=150)

    # --- Filter low-coverage reads (typical first step on real data) ---
    prm_filtered = prm.filter_reads(min_sites=8)
    keep = np.isin(prm.read_ids, prm_filtered.read_ids)
    true_labels = true_labels[keep]
    fragment_lengths = fragment_lengths[keep]
    print(f"  {prm_filtered.n_reads} reads remain after coverage filter")

    # --- Classify ---
    print("[example] Classifying single-molecule footprint states")
    result = classify.classify_footprints(prm_filtered, feature_center=cfg.feature_center)
    print("  state counts:", result.counts())

    # Confusion against ground truth.
    df = pd.DataFrame({"true": true_labels, "predicted": result.states})
    confusion = pd.crosstab(df["true"], df["predicted"], dropna=False)
    print("[example] Confusion matrix (rows = truth, cols = prediction):")
    print(confusion)

    # --- Sort by motif accessibility for the heatmap ---
    sorted_prm = per_read.sort_by_pattern(prm_filtered, feature_center=cfg.feature_center, window=80)

    # --- Plots ---
    print("[example] Rendering plots")
    fig_hm = viz.single_molecule_heatmap(
        sorted_prm,
        feature_center=cfg.feature_center,
        title="Synthetic SMF -- single-molecule heatmap (sorted by motif accessibility)",
    )
    fig_hm.savefig(out_dir / "single_molecule_heatmap.png", dpi=150)

    fig_avg = viz.average_accessibility_plot(
        prm_filtered,
        feature_center=cfg.feature_center,
        smooth_bp=20,
        title="Synthetic SMF -- average GpC accessibility",
    )
    fig_avg.savefig(out_dir / "average_accessibility.png", dpi=150)

    fig_states = viz.state_composition_plot(result.counts())
    fig_states.savefig(out_dir / "state_composition.png", dpi=150)

    fig_cov = viz.coverage_per_position_plot(
        prm_filtered, feature_center=cfg.feature_center, low_coverage_threshold=20,
        title="Synthetic SMF -- coverage per informative GpC site",
    )
    fig_cov.savefig(out_dir / "coverage_per_position.png", dpi=150)

    # --- Co-occupancy: motif vs. immediate flank ---
    print("[example] Co-occupancy (motif vs. +60bp downstream flank)")
    co = stats.co_occupancy(
        prm_filtered,
        bin_a=(cfg.feature_center - 15, cfg.feature_center + 15),
        bin_b=(cfg.feature_center + 45, cfg.feature_center + 75),
    )
    print(f"  n={co.n}, p11={co.p11:.3f}, p10={co.p10:.3f}, "
          f"p01={co.p01:.3f}, p00={co.p00:.3f}, log2_OR={co.log2_or:.2f}")

    # --- Persist a summary table for downstream / re-analysis ---
    summary = pd.DataFrame({
        "read_id": prm_filtered.read_ids,
        "true_state": true_labels,
        "predicted_state": result.states,
        "motif_score": result.motif_score,
        "inner_flank_score": result.inner_flank_score,
        "outer_flank_score": result.outer_flank_score,
    })
    summary.to_csv(out_dir / "per_read_summary.csv", index=False)
    confusion.to_csv(out_dir / "confusion_matrix.csv")

    print(f"[example] Outputs written to {out_dir}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
