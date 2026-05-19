"""Plotting helpers for SMF analysis.

All functions return the matplotlib ``Figure`` so they're easy to compose, save, and
test. They use matplotlib only -- no seaborn, no extra style files.
"""
from __future__ import annotations

from typing import Sequence

import matplotlib
import numpy as np

# Use a non-interactive backend so this module imports cleanly in headless / CI envs.
if matplotlib.get_backend().lower() != "agg":
    matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.figure import Figure

from .per_read import PerReadMatrix


def single_molecule_heatmap(
    prm: PerReadMatrix,
    *,
    feature_center: int | None = None,
    title: str = "Single-molecule SMF heatmap",
    cmap_methylated: str = "tab:blue",     # accessible / methylated GpC
    cmap_unmethylated: str = "tab:red",    # protected / unmethylated GpC -- field convention
    nan_color: str = "white",
    figsize: tuple[float, float] = (8.0, 6.0),
) -> Figure:
    """Render the canonical SMF single-molecule heatmap.

    Each row is a read, each column an informative GpC site. Methylated calls
    (accessible) are dark, unmethylated calls (protected) are light, missing calls
    are white. A vertical line marks ``feature_center`` if given.
    """
    cmap = ListedColormap([nan_color, cmap_unmethylated, cmap_methylated])
    norm = BoundaryNorm([-1.5, -0.5, 0.5, 1.5], cmap.N)

    data = np.where(np.isnan(prm.matrix), -1.0, prm.matrix)

    fig, ax = plt.subplots(figsize=figsize)
    extent = (
        float(prm.sites.min()),
        float(prm.sites.max()),
        prm.n_reads,
        0,
    )
    ax.imshow(data, aspect="auto", interpolation="nearest", cmap=cmap, norm=norm, extent=extent)
    ax.set_xlabel("Reference position (bp)")
    ax.set_ylabel(f"Reads (n={prm.n_reads})")
    ax.set_title(title)

    if feature_center is not None:
        ax.axvline(feature_center, color="black", linewidth=0.8, linestyle="--")

    # Build a small custom legend rather than a colourbar -- it's three categories.
    handles = [
        plt.Rectangle((0, 0), 1, 1, color=cmap_methylated, label="Accessible (methylated GpC)"),
        plt.Rectangle((0, 0), 1, 1, color=cmap_unmethylated, label="Protected (unmethylated GpC)"),
        plt.Rectangle((0, 0), 1, 1, color=nan_color, ec="grey", label="No call"),
    ]
    ax.legend(handles=handles, loc="upper right", fontsize=8, frameon=True)
    fig.tight_layout()
    return fig


def average_accessibility_plot(
    prm: PerReadMatrix,
    *,
    feature_center: int | None = None,
    smooth_bp: int = 20,
    title: str = "Average accessibility",
    figsize: tuple[float, float] = (8.0, 3.0),
) -> Figure:
    """Plot the smoothed average GpC methylation across the window."""
    from .stats import average_accessibility_profile

    grid, smoothed = average_accessibility_profile(prm, smooth_bp=smooth_bp)

    fig, ax = plt.subplots(figsize=figsize)
    ax.plot(grid, smoothed, color="tab:blue", linewidth=1.4)
    ax.set_ylim(0, 1)
    ax.set_xlabel("Reference position (bp)")
    ax.set_ylabel("Mean GpC methylation")
    ax.set_title(title)
    if feature_center is not None:
        ax.axvline(feature_center, color="black", linewidth=0.8, linestyle="--")
    fig.tight_layout()
    return fig


def single_molecule_heatmap_filled(
    prm: PerReadMatrix,
    *,
    feature_center: int | None = None,
    title: str = "Single-molecule SMF heatmap (filled)",
    color_protected: str = "tab:red",
    color_accessible: str = "tab:blue",
    nan_color: str = "white",
    figsize: tuple[float, float] = (8.0, 6.0),
) -> Figure:
    """Continuous-style SMF heatmap matching the field's convention.

    Unlike :func:`single_molecule_heatmap`, this fills the space *between* informative
    GpC sites. For each read, every base in the per-base grid is colored by the
    methylation status of the **nearest informative call**, bounded by the read's
    leftmost and rightmost calls. Outside that range the cell is left white.

    Color convention (matches the standard SMF literature):

    * **Red**   = unmethylated GpC  -> protected DNA (e.g. wrapped around a nucleosome)
    * **Blue**  = methylated GpC    -> accessible DNA (linker / open chromatin)
    * **White** = the read had no informative calls in this region

    The continuous fill makes nucleosome and TF footprints visually obvious as solid
    red blocks; isolated red dots in a sea of blue look like noise.
    """
    if prm.n_reads == 0 or prm.n_sites == 0:
        raise ValueError("empty per-read matrix; nothing to plot")

    pos_min = int(prm.sites.min())
    pos_max = int(prm.sites.max())
    grid = np.arange(pos_min, pos_max + 1)
    filled = np.full((prm.n_reads, len(grid)), np.nan, dtype=float)

    for r_idx in range(prm.n_reads):
        row = prm.matrix[r_idx, :]
        valid = ~np.isnan(row)
        if not valid.any():
            continue
        valid_sites = prm.sites[valid]
        valid_vals = row[valid]

        # For each base in the grid, find the index of the nearest informative call.
        idx_right = np.searchsorted(valid_sites, grid)
        idx_left = np.clip(idx_right - 1, 0, len(valid_sites) - 1)
        idx_right = np.clip(idx_right, 0, len(valid_sites) - 1)
        dist_left = np.abs(grid - valid_sites[idx_left])
        dist_right = np.abs(grid - valid_sites[idx_right])
        nearest = np.where(dist_left <= dist_right, idx_left, idx_right)
        row_filled = valid_vals[nearest]

        # Bound by the read's actual coverage so we don't extrapolate beyond it.
        in_range = (grid >= valid_sites[0]) & (grid <= valid_sites[-1])
        filled[r_idx] = np.where(in_range, row_filled, np.nan)

    # Render using a 3-bin colormap: NaN -> white, 0 -> red, 1 -> blue.
    data = np.where(np.isnan(filled), -1.0, filled)
    cmap = ListedColormap([nan_color, color_protected, color_accessible])
    norm = BoundaryNorm([-1.5, -0.5, 0.5, 1.5], cmap.N)

    fig, ax = plt.subplots(figsize=figsize)
    extent = (float(pos_min), float(pos_max), prm.n_reads, 0)
    ax.imshow(data, aspect="auto", interpolation="nearest", cmap=cmap, norm=norm, extent=extent)
    ax.set_xlabel("Reference position (bp)")
    ax.set_ylabel(f"Reads (n={prm.n_reads})")
    ax.set_title(title)
    if feature_center is not None:
        ax.axvline(feature_center, color="black", linewidth=0.8, linestyle="--")

    handles = [
        plt.Rectangle((0, 0), 1, 1, color=color_accessible, label="Accessible (methylated GpC)"),
        plt.Rectangle((0, 0), 1, 1, color=color_protected, label="Protected (unmethylated GpC)"),
        plt.Rectangle((0, 0), 1, 1, color=nan_color, ec="grey", label="No call"),
    ]
    ax.legend(handles=handles, loc="upper right", fontsize=8, frameon=True)
    fig.tight_layout()
    return fig


def coverage_per_position_plot(
    prm: PerReadMatrix,
    *,
    feature_center: int | None = None,
    title: str = "Sequencing coverage per informative cytosine",
    figsize: tuple[float, float] = (8.0, 3.0),
    low_coverage_threshold: int | None = None,
) -> Figure:
    """Bar plot of how many reads have an informative call at each cytosine.

    Useful for spotting positions you should *not* trust in the heatmap because only a
    handful of molecules cover them. If ``low_coverage_threshold`` is given, sites
    below that threshold are highlighted in red.
    """
    cov = prm.coverage_per_site()
    fig, ax = plt.subplots(figsize=figsize)

    colors = ["tab:blue"] * len(cov)
    if low_coverage_threshold is not None:
        for i, c in enumerate(cov):
            if c < low_coverage_threshold:
                colors[i] = "tab:red"

    ax.bar(prm.sites, cov, color=colors, width=max(1, int((prm.sites.max() - prm.sites.min()) / max(1, len(prm.sites)))))
    ax.set_xlabel("Reference position (bp)")
    ax.set_ylabel("Reads covering site")
    ax.set_title(title)
    if feature_center is not None:
        ax.axvline(feature_center, color="black", linewidth=0.8, linestyle="--")
    if low_coverage_threshold is not None:
        ax.axhline(low_coverage_threshold, color="grey", linewidth=0.6, linestyle=":")
    fig.tight_layout()
    return fig


def fragment_length_plot(
    lengths: np.ndarray | Sequence[int],
    *,
    bins: int = 60,
    min_keep: int | None = None,
    title: str = "Fragment length distribution",
    figsize: tuple[float, float] = (6.0, 3.5),
) -> Figure:
    """Histogram of merged-read / fragment lengths.

    For SMF this is a critical QC plot -- short fragments usually mean bisulfite damage
    or incomplete amplicon, and they distort the heatmap. ``min_keep`` draws a vertical
    line at your length filter (e.g. 500 bp) so you can see how much you're discarding.
    """
    arr = np.asarray(lengths, dtype=int)
    fig, ax = plt.subplots(figsize=figsize)
    ax.hist(arr, bins=bins, color="tab:blue", edgecolor="white")
    ax.set_xlabel("Read / fragment length (bp)")
    ax.set_ylabel("Count")
    ax.set_title(title)
    median = int(np.median(arr)) if len(arr) else 0
    ax.axvline(median, color="black", linewidth=0.8, linestyle="--", label=f"median = {median}")
    if min_keep is not None:
        n_drop = int(np.sum(arr < min_keep))
        ax.axvline(min_keep, color="tab:red", linewidth=0.8, linestyle="--",
                   label=f"min_keep = {min_keep} ({n_drop} below)")
    ax.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    return fig


def state_composition_plot(
    counts: dict[str, int],
    *,
    order: Sequence[str] = ("accessible", "tf_bound", "nucleosome", "unclassified"),
    colors: dict[str, str] | None = None,
    title: str = "Single-molecule footprint state composition",
    figsize: tuple[float, float] = (5.0, 3.0),
) -> Figure:
    """Stacked-bar / pie summary of footprint state fractions.

    Renders a single horizontal stacked bar -- much more readable than a pie when there
    are 3-4 categories.
    """
    palette = {
        "accessible": "tab:green",
        "tf_bound": "tab:red",
        "nucleosome": "tab:purple",
        "unclassified": "lightgrey",
    }
    if colors:
        palette.update(colors)

    total = max(sum(counts.values()), 1)
    fig, ax = plt.subplots(figsize=figsize)

    left = 0.0
    for state in order:
        n = counts.get(state, 0)
        if n == 0:
            continue
        frac = n / total
        ax.barh(0, frac, left=left, color=palette.get(state, "tab:grey"), edgecolor="white", label=f"{state} (n={n})")
        if frac > 0.04:
            ax.text(left + frac / 2, 0, f"{frac*100:.0f}%", ha="center", va="center", fontsize=9, color="white")
        left += frac

    ax.set_xlim(0, 1)
    ax.set_yticks([])
    ax.set_xlabel("Fraction of molecules")
    ax.set_title(title)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.4), ncol=2, fontsize=8, frameon=False)
    fig.tight_layout()
    return fig
