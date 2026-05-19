"""Single-molecule footprint state classification.

Given a per-read GpC (accessibility) matrix centred on a feature (e.g. a TF motif),
classify each molecule into one of three biologically meaningful states:

* ``ACCESSIBLE`` -- high accessibility throughout the window
* ``TF_BOUND``   -- a small (~30 bp) protected patch over the motif, surrounded by
                     accessible flanks
* ``NUCLEOSOME`` -- a wide (~150 bp) protected patch covering the motif and most of
                     the surrounding window

This is the "three-state" classifier from Krebs/Sönmezer-style analyses, simplified to
one short read pair per molecule. It is intentionally rule-based and inspectable;
swap in a probabilistic model if your data warrants it.
"""
from __future__ import annotations

import warnings
from dataclasses import dataclass
from enum import Enum

import numpy as np

from .per_read import PerReadMatrix


class FootprintState(str, Enum):
    ACCESSIBLE = "accessible"
    TF_BOUND = "tf_bound"
    NUCLEOSOME = "nucleosome"
    UNCLASSIFIED = "unclassified"


@dataclass
class ClassifierConfig:
    """Tunable thresholds for the three-state classifier.

    All windows are *symmetric* and given as half-widths in bp.

    Defaults are set for a vertebrate TF-binding-site analysis with M.CviPI GpC data;
    inspect your average accessibility plot first and adjust if your numbers look off.
    """

    motif_half_width: int = 15        # +/- around feature_center -> "motif" core
    flank_inner: int = 35             # core flank for TF-bound test (just outside motif)
    flank_outer: int = 100            # broader flank for nucleosome test
    min_motif_calls: int = 2          # require at least N informative sites in the core
    min_flank_calls: int = 2          # required flank sites for either flank
    motif_protected_max: float = 0.33 # mean accessibility in motif <= this => "protected"
    flank_accessible_min: float = 0.6 # mean accessibility in inner flank >= this => "accessible"
    nucleosome_protected_max: float = 0.4  # mean accessibility across outer flank for nuc.


@dataclass
class ClassificationResult:
    states: np.ndarray              # (n_reads,) of FootprintState string values
    motif_score: np.ndarray         # (n_reads,) mean meth in motif (NaN -> insufficient calls)
    inner_flank_score: np.ndarray   # (n_reads,) mean meth in inner flank
    outer_flank_score: np.ndarray   # (n_reads,) mean meth in outer flank

    def counts(self) -> dict[str, int]:
        unique, counts = np.unique(self.states, return_counts=True)
        return {str(u): int(c) for u, c in zip(unique, counts)}

    def fractions(self) -> dict[str, float]:
        n = max(len(self.states), 1)
        return {k: v / n for k, v in self.counts().items()}


def _window_mean(
    matrix: np.ndarray,
    sites: np.ndarray,
    center: int,
    half_width: int,
    min_calls: int,
) -> np.ndarray:
    """Per-read mean methylation in ``[center-half_width, center+half_width]``,
    NaN where a read has fewer than ``min_calls`` informative sites in the window.
    """
    mask = (sites >= center - half_width) & (sites <= center + half_width)
    if not mask.any():
        return np.full(matrix.shape[0], np.nan)
    sub = matrix[:, mask]
    n_calls = np.sum(~np.isnan(sub), axis=1)
    with np.errstate(invalid="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        means = np.nanmean(sub, axis=1)
    means = np.where(n_calls >= min_calls, means, np.nan)
    return means


def classify_footprints(
    prm: PerReadMatrix,
    feature_center: int,
    config: ClassifierConfig | None = None,
) -> ClassificationResult:
    """Classify each read in ``prm`` into a :class:`FootprintState` around
    ``feature_center`` (a 0-based reference coordinate)."""
    cfg = config or ClassifierConfig()

    motif = _window_mean(prm.matrix, prm.sites, feature_center, cfg.motif_half_width, cfg.min_motif_calls)
    inner = _outer_minus_inner_mean(
        prm.matrix, prm.sites, feature_center,
        inner_half=cfg.motif_half_width, outer_half=cfg.flank_inner,
        min_calls=cfg.min_flank_calls,
    )
    outer = _outer_minus_inner_mean(
        prm.matrix, prm.sites, feature_center,
        inner_half=cfg.motif_half_width, outer_half=cfg.flank_outer,
        min_calls=cfg.min_flank_calls,
    )

    states = np.full(prm.n_reads, FootprintState.UNCLASSIFIED.value, dtype=object)

    # Order matters: classify the most specific patterns first.
    tf_bound = (
        ~np.isnan(motif) & ~np.isnan(inner)
        & (motif <= cfg.motif_protected_max)
        & (inner >= cfg.flank_accessible_min)
    )
    nucleosome = (
        ~np.isnan(motif) & ~np.isnan(outer)
        & (motif <= cfg.motif_protected_max)
        & (outer <= cfg.nucleosome_protected_max)
        & ~tf_bound
    )
    accessible = (
        ~np.isnan(motif) & ~np.isnan(outer)
        & (motif > cfg.motif_protected_max)
        & (outer > cfg.nucleosome_protected_max)
    )

    states[tf_bound] = FootprintState.TF_BOUND.value
    states[nucleosome] = FootprintState.NUCLEOSOME.value
    states[accessible] = FootprintState.ACCESSIBLE.value

    return ClassificationResult(
        states=states,
        motif_score=motif,
        inner_flank_score=inner,
        outer_flank_score=outer,
    )


def _outer_minus_inner_mean(
    matrix: np.ndarray,
    sites: np.ndarray,
    center: int,
    *,
    inner_half: int,
    outer_half: int,
    min_calls: int,
) -> np.ndarray:
    """Mean methylation over the **annulus** ``[-outer_half, -inner_half] U [inner_half, outer_half]``
    around ``center``, per read. NaN if not enough informative sites."""
    mask_outer = (sites >= center - outer_half) & (sites <= center + outer_half)
    mask_inner = (sites >= center - inner_half) & (sites <= center + inner_half)
    annulus = mask_outer & ~mask_inner
    if not annulus.any():
        return np.full(matrix.shape[0], np.nan)
    sub = matrix[:, annulus]
    n_calls = np.sum(~np.isnan(sub), axis=1)
    with np.errstate(invalid="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        means = np.nanmean(sub, axis=1)
    return np.where(n_calls >= min_calls, means, np.nan)
