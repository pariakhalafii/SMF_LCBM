"""Statistical summaries over per-read methylation matrices."""
from __future__ import annotations

import warnings
from dataclasses import dataclass

import numpy as np

from .per_read import PerReadMatrix


@dataclass
class CoOccupancy:
    """Pairwise co-methylation between two sets of sites.

    ``p11`` is the fraction of reads (with calls in both bins) that are methylated at
    both bins; ``p10``/``p01``/``p00`` are the other three combinations.

    The log odds ratio :math:`\\log_2 \\frac{p_{11} p_{00}}{p_{10} p_{01}}` is a
    standard summary of co-occupancy on the same molecule.
    """
    n: int
    p11: float
    p10: float
    p01: float
    p00: float

    @property
    def log2_or(self) -> float:
        num = max(self.p11 * self.p00, 1e-9)
        den = max(self.p10 * self.p01, 1e-9)
        return float(np.log2(num / den))


def _bin_mean(matrix: np.ndarray, mask: np.ndarray, threshold: float = 0.5) -> np.ndarray:
    """Per-read average methylation over the columns in ``mask``, then thresholded
    to a 0/1 call. Returns NaN where a read has no calls in the bin."""
    if not mask.any():
        return np.full(matrix.shape[0], np.nan)
    sub = matrix[:, mask]
    has_call = np.any(~np.isnan(sub), axis=1)
    with np.errstate(invalid="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        means = np.nanmean(sub, axis=1)
    binarised = np.where(means >= threshold, 1.0, 0.0)
    return np.where(has_call, binarised, np.nan)


def co_occupancy(
    prm: PerReadMatrix,
    bin_a: tuple[int, int],
    bin_b: tuple[int, int],
    threshold: float = 0.5,
) -> CoOccupancy:
    """Co-methylation summary between two reference-coordinate intervals.

    ``bin_a`` and ``bin_b`` are ``(start, end)`` intervals in the same coordinate
    system as ``prm.sites``. A read contributes only if it has at least one informative
    call in both bins.
    """
    mask_a = (prm.sites >= bin_a[0]) & (prm.sites <= bin_a[1])
    mask_b = (prm.sites >= bin_b[0]) & (prm.sites <= bin_b[1])

    a = _bin_mean(prm.matrix, mask_a, threshold)
    b = _bin_mean(prm.matrix, mask_b, threshold)
    keep = ~np.isnan(a) & ~np.isnan(b)
    a, b = a[keep], b[keep]

    n = int(keep.sum())
    if n == 0:
        return CoOccupancy(n=0, p11=float("nan"), p10=float("nan"), p01=float("nan"), p00=float("nan"))

    p11 = float(np.mean((a == 1) & (b == 1)))
    p10 = float(np.mean((a == 1) & (b == 0)))
    p01 = float(np.mean((a == 0) & (b == 1)))
    p00 = float(np.mean((a == 0) & (b == 0)))
    return CoOccupancy(n=n, p11=p11, p10=p10, p01=p01, p00=p00)


def average_accessibility_profile(
    prm: PerReadMatrix,
    *,
    smooth_bp: int | None = 20,
) -> tuple[np.ndarray, np.ndarray]:
    """Return ``(positions, mean_methylation)`` ready for plotting.

    If ``smooth_bp`` is given, a simple sliding-window mean is applied over base-pair
    space (not site space), which produces a much more readable curve when site density
    is uneven.
    """
    sites = prm.sites
    means = prm.average_accessibility()
    if smooth_bp is None or len(sites) < 3:
        return sites, means

    # Project onto a regular bp grid before smoothing so unevenly spaced sites don't
    # bias the moving average.
    grid = np.arange(int(sites.min()), int(sites.max()) + 1)
    interp = np.interp(grid, sites, np.where(np.isnan(means), 0.5, means))
    half = max(1, smooth_bp // 2)
    kernel = np.ones(2 * half + 1) / (2 * half + 1)
    # Pad with edge values (not zeros) so the smoothing doesn't pull edge positions
    # toward 0.5 and produce a false SMF spike at the start/end of the window.
    padded = np.pad(interp, half, mode="edge")
    smoothed = np.convolve(padded, kernel, mode="valid")
    smoothed = np.clip(smoothed, 0.0, 1.0)
    return grid, smoothed
