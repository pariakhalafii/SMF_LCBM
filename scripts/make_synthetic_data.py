"""Generate a synthetic per-read GpC methylation matrix for SMF testing.

Models three molecule populations across a window centred on a TF motif:

* fully accessible -- high GpC methylation everywhere
* TF-bound        -- ~30 bp protected patch over the motif, accessible elsewhere
* nucleosome      -- ~150 bp protected patch covering the motif

The output is a :class:`smf.per_read.PerReadMatrix` ready for downstream analysis.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from smf.per_read import PerReadMatrix


@dataclass
class SyntheticConfig:
    n_reads: int = 600
    window_bp: int = 600                     # window is [-window_bp/2, +window_bp/2] around centre
    feature_center: int = 1000               # arbitrary reference coord
    gpc_density_bp: int = 8                  # average bp between informative GpC sites
    fraction_accessible: float = 0.45
    fraction_tf_bound: float = 0.20
    fraction_nucleosome: float = 0.30        # remainder is "ambiguous / partial"
    # TF footprint half-width. Realistic TF footprints are 30-50 bp wide TOTAL, so a
    # half-width around 15 bp lets the classifier's inner-flank annulus pick up
    # accessible flanks immediately outside the motif core.
    tf_footprint_bp: int = 15
    nucleosome_footprint_bp: int = 75        # half-width of nucleosome protection
    p_meth_accessible: float = 0.85          # P(methylation | accessible) for M.CviPI
    p_meth_protected: float = 0.10           # P(methylation | protected) -- background
    dropout_rate: float = 0.10               # P(no call) per read x site
    read_span_bp: int = 250                  # reads cover only a chunk of the window
    seed: int = 0


def make_synthetic(
    cfg: SyntheticConfig | None = None,
) -> tuple[PerReadMatrix, np.ndarray, np.ndarray]:
    """Return ``(per_read_matrix, true_labels, fragment_lengths)``.

    ``true_labels`` is a string array of length ``n_reads`` giving the ground-truth
    state for each read, useful for evaluating the classifier.
    ``fragment_lengths`` is the simulated merged-read length per molecule (bp), for
    fragment-length QC plots.
    """
    cfg = cfg or SyntheticConfig()
    rng = np.random.default_rng(cfg.seed)

    half = cfg.window_bp // 2
    win_start = cfg.feature_center - half
    win_end = cfg.feature_center + half

    # Place GpC sites with mild jitter so density is roughly cfg.gpc_density_bp.
    n_sites_target = max(2, cfg.window_bp // cfg.gpc_density_bp)
    base_positions = np.linspace(win_start + 2, win_end - 2, n_sites_target)
    jitter = rng.integers(-cfg.gpc_density_bp // 2, cfg.gpc_density_bp // 2, size=n_sites_target)
    sites = np.unique(np.clip(base_positions.astype(int) + jitter, win_start, win_end))

    # Decide each read's true state.
    fracs = np.array([cfg.fraction_accessible, cfg.fraction_tf_bound, cfg.fraction_nucleosome])
    fracs = np.append(fracs, max(0.0, 1.0 - fracs.sum()))
    fracs = fracs / fracs.sum()
    labels = rng.choice(
        ["accessible", "tf_bound", "nucleosome", "unclassified"],
        size=cfg.n_reads,
        p=fracs,
    )

    # Per-site protection probability per read.
    matrix = np.full((cfg.n_reads, len(sites)), np.nan, dtype=float)
    fragment_lengths = np.zeros(cfg.n_reads, dtype=int)
    for r in range(cfg.n_reads):
        # Pick a random read interval inside the window so coverage is realistic.
        # Add a little length jitter (~+/-15%) to mimic real fragment-length spread.
        nominal = min(cfg.read_span_bp, cfg.window_bp)
        span = max(40, int(nominal * rng.uniform(0.6, 1.0)))
        fragment_lengths[r] = span
        read_start = rng.integers(win_start, win_end - span + 1)
        read_end = read_start + span
        covered = (sites >= read_start) & (sites <= read_end)

        protected = _protected_mask(sites, labels[r], cfg)
        p_meth = np.where(protected, cfg.p_meth_protected, cfg.p_meth_accessible)

        for j in np.where(covered)[0]:
            if rng.random() < cfg.dropout_rate:
                continue   # no call -> NaN
            matrix[r, j] = float(rng.random() < p_meth[j])

    prm = PerReadMatrix(
        matrix=matrix,
        sites=sites,
        read_ids=[f"read_{i:05d}" for i in range(cfg.n_reads)],
        chrom="synth",
        region_start=win_start,
        region_end=win_end,
    )
    return prm, labels, fragment_lengths


def _protected_mask(sites: np.ndarray, label: str, cfg: SyntheticConfig) -> np.ndarray:
    """Boolean mask of which sites are inside the protection footprint for ``label``."""
    c = cfg.feature_center
    if label == "accessible":
        return np.zeros_like(sites, dtype=bool)
    if label == "tf_bound":
        return (sites >= c - cfg.tf_footprint_bp) & (sites <= c + cfg.tf_footprint_bp)
    if label == "nucleosome":
        return (sites >= c - cfg.nucleosome_footprint_bp) & (sites <= c + cfg.nucleosome_footprint_bp)
    # "unclassified" -- random scattered protection
    return np.zeros_like(sites, dtype=bool)
