"""smf -- Single Molecule Footprinting analysis utilities (bisulfite / Illumina).

The package is organised as small, composable modules:

* :mod:`smf.context`  -- determine sequence context (CpG / GpC / GCG) for a position
* :mod:`smf.io`       -- read Bismark CX reports and per-read methylation from BAMs
* :mod:`smf.per_read` -- build per-read methylation matrices around regions of interest
* :mod:`smf.classify` -- classify single molecules into footprint states
* :mod:`smf.stats`    -- co-occupancy, average accessibility, summary stats
* :mod:`smf.viz`      -- single-molecule heatmaps and average plots

The downstream modules (``per_read``, ``classify``, ``stats``, ``viz``) operate on a
plain ``numpy`` per-read methylation matrix and don't require any external bioinformatics
tools, which makes them easy to test and reuse.
"""

from . import classify, context, io, per_read, stats, viz  # noqa: F401

__all__ = ["context", "io", "per_read", "classify", "stats", "viz"]

__version__ = "0.1.0"
