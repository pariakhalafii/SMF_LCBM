"""Per-read methylation matrix utilities.

The central data structure for downstream SMF analysis is a per-read methylation matrix:

* rows  = molecules (reads / read-pairs)
* cols  = informative cytosines in the region of interest, ordered by position
* values = ``1.0`` (methylated), ``0.0`` (unmethylated), or ``np.nan`` (no call).

For accessibility analysis only **GCH** sites should be in the matrix; for endogenous
methylation use **HCG**. Mixing them is almost always a bug -- keep two matrices.
"""
from __future__ import annotations

import re
import warnings
from dataclasses import dataclass
from typing import Sequence

import numpy as np


@dataclass
class PerReadMatrix:
    """A per-read methylation matrix bundled with its column coordinates."""

    matrix: np.ndarray         # (n_reads, n_sites), values in {0, 1, NaN}
    sites: np.ndarray          # (n_sites,) reference positions, 0-based
    read_ids: list[str]
    chrom: str = ""
    region_start: int = 0      # 0-based inclusive
    region_end: int = 0        # 0-based exclusive

    def __post_init__(self) -> None:
        if self.matrix.shape != (len(self.read_ids), len(self.sites)):
            raise ValueError(
                f"matrix shape {self.matrix.shape} doesn't match "
                f"({len(self.read_ids)}, {len(self.sites)})"
            )

    @property
    def n_reads(self) -> int:
        return self.matrix.shape[0]

    @property
    def n_sites(self) -> int:
        return self.matrix.shape[1]

    def coverage_per_read(self) -> np.ndarray:
        """Number of informative (non-NaN) sites per read."""
        return np.sum(~np.isnan(self.matrix), axis=1)

    def coverage_per_site(self) -> np.ndarray:
        """Number of reads that cover each site."""
        return np.sum(~np.isnan(self.matrix), axis=0)

    def filter_reads(self, min_sites: int) -> "PerReadMatrix":
        """Return a new matrix keeping only reads with ``>= min_sites`` informative
        calls. Useful before clustering or classification."""
        keep = self.coverage_per_read() >= min_sites
        return PerReadMatrix(
            matrix=self.matrix[keep],
            sites=self.sites,
            read_ids=[rid for rid, k in zip(self.read_ids, keep) if k],
            chrom=self.chrom,
            region_start=self.region_start,
            region_end=self.region_end,
        )

    def average_accessibility(self) -> np.ndarray:
        """Per-site mean methylation, ignoring NaNs. NaN where no read covers a site."""
        with np.errstate(invalid="ignore"):
            return np.nanmean(self.matrix, axis=0)


def per_read_matrix_from_calls(
    calls,                         # Iterable[PerReadCall]
    sites: Sequence[int],
    *,
    chrom: str = "",
    region_start: int = 0,
    region_end: int = 0,
) -> PerReadMatrix:
    """Convenience wrapper: build a :class:`PerReadMatrix` directly from
    :class:`smf.io.PerReadCall` objects."""
    from .io import calls_to_matrix

    matrix, read_ids = calls_to_matrix(list(calls), list(sites))
    return PerReadMatrix(
        matrix=matrix,
        sites=np.asarray(sites, dtype=int),
        read_ids=read_ids,
        chrom=chrom,
        region_start=region_start,
        region_end=region_end,
    )


def sort_by_pattern(
    prm: PerReadMatrix,
    *,
    feature_center: int | None = None,
    window: int = 100,
) -> PerReadMatrix:
    """Sort reads so that the most-protected molecules at the feature centre cluster
    together at the top -- this is the canonical layout for SMF heatmaps.

    Strategy: rank each read by its mean accessibility within ``[center-window, center+window]``,
    ascending (most protected first). Reads with no calls in the window go to the bottom.
    """
    if feature_center is None:
        feature_center = int((prm.region_start + prm.region_end) / 2)
    in_window = (prm.sites >= feature_center - window) & (prm.sites <= feature_center + window)
    if not in_window.any():
        return prm
    sub = prm.matrix[:, in_window]
    with np.errstate(invalid="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        score = np.nanmean(sub, axis=1)
    score_filled = np.where(np.isnan(score), np.inf, score)
    order = np.argsort(score_filled, kind="stable")
    return PerReadMatrix(
        matrix=prm.matrix[order],
        sites=prm.sites,
        read_ids=[prm.read_ids[i] for i in order],
        chrom=prm.chrom,
        region_start=prm.region_start,
        region_end=prm.region_end,
    )


# ---------------------------------------------------------------------------
# Pair-merging for paired-end Bismark BAMs
# ---------------------------------------------------------------------------

# Trailing Illumina-style R1/R2 indicators that distinguish the two mates of a
# pair within a BAM. We strip these to obtain a per-MOLECULE identifier:
#   * `/1` / `/2`                      (older SAM convention)
#   * `_1:N:0:CGAACAAG`                (Illumina comment glued in by some tools)
#   * ` 1:N:0:CGAACAAG`                (with space preserved)
#   * `_1` / `_2` at the very end       (minimal form)
_PAIR_SUFFIX_RE = re.compile(r"[_/ ][12](?::N:\d+:[A-Za-z0-9+]+)?$")


def _strip_pair_suffix(read_id: str) -> str:
    """Strip the trailing R1/R2 indicator from a Bismark paired-end read name."""
    return _PAIR_SUFFIX_RE.sub("", read_id)


def merge_paired_reads(prm: "PerReadMatrix") -> "PerReadMatrix":
    """Merge R1 and R2 mates of paired-end reads into one row per DNA molecule.

    In a paired-end Bismark BAM, R1 and R2 of the same molecule appear as two
    separate records and sometimes end up with slightly different ``read_id``
    strings (Illumina mate indicators get attached to the name). This function
    strips those indicators and combines rows whose normalized read_id matches,
    yielding one row per physical DNA molecule.

    Combining rule per position:
        * if only one mate has a call -> use it
        * if both mates have a call and they agree -> use that value
        * if both have a call and disagree -> 0.5 (rendered as methylated by the
          heatmap's BoundaryNorm; conflicts are rare and indicate sequencing noise
          in the overlap region between R1 and R2)

    Returns a new :class:`PerReadMatrix` with one row per molecule.
    """
    import pandas as pd  # local import to keep base module light

    if prm.n_reads == 0:
        return prm

    # Group rows by their normalized (R1/R2-stripped) identifier.
    norm_ids = [_strip_pair_suffix(rid) for rid in prm.read_ids]
    df = pd.DataFrame(prm.matrix, index=norm_ids, columns=prm.sites)
    # groupby.mean() ignores NaNs by default (skipna=True). So a position covered
    # only by one mate will keep that mate's value, and a position covered by
    # neither stays NaN.
    merged = df.groupby(level=0, sort=False).mean()

    return PerReadMatrix(
        matrix=merged.to_numpy(dtype=float),
        sites=prm.sites,
        read_ids=list(merged.index),
        chrom=prm.chrom,
        region_start=prm.region_start,
        region_end=prm.region_end,
    )
