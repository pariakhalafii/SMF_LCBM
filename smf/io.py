"""I/O for Bismark CX reports and per-read methylation calls.

We support two common SMF inputs:

1. **CX report** (``bismark_methylation_extractor --CX --cytosine_report``):
   a TSV of ``chrom  pos  strand  count_meth  count_unmeth  context  trinucleotide``.
   Useful for *bulk* per-cytosine methylation; loaded with :func:`read_cx_report`.

2. **Per-read methylation from a Bismark BAM** via the ``XM`` tag:
   each aligned read carries a per-base call string where ``Z``/``z`` = methylated /
   unmethylated CpG, ``X``/``x`` = CHG, ``H``/``h`` = CHH, and ``U``/``u`` = unknown.
   For NOMe-seq data with ``--CX``, GpC calls live in the same XM string but with the
   correct context resolved from the reference. We re-resolve context from the
   reference here to be safe (Bismark's older versions are inconsistent across CX/non-CX
   modes).

Reading a real BAM requires :mod:`pysam`; reading the FASTA requires :mod:`pyfaidx`.
Both are optional dependencies -- the rest of the package works on plain numpy
matrices, which is how the synthetic example exercises everything end-to-end.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, Sequence

import numpy as np
import pandas as pd

from .context import Context, context_at


# ---------------------------------------------------------------------------
# CX report (bulk)
# ---------------------------------------------------------------------------

CX_COLUMNS = [
    "chrom",
    "pos",          # 1-based in Bismark output
    "strand",
    "count_meth",
    "count_unmeth",
    "context",      # Bismark labels: CpG / CHG / CHH (does NOT split GpC out)
    "trinucleotide",
]


def read_cx_report(path: str | Path) -> pd.DataFrame:
    """Read a Bismark CX report (gzipped or plain) into a DataFrame.

    Adds a ``smf_context`` column with the resolved :class:`smf.context.Context` for
    each cytosine, and a ``frac_meth`` column = ``count_meth / (count_meth+count_unmeth)``
    (NaN where coverage is 0).
    """
    path = Path(path)
    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=CX_COLUMNS,
        dtype={"chrom": str, "strand": str, "context": str, "trinucleotide": str},
        compression="infer",
    )

    df["smf_context"] = df["trinucleotide"].apply(_resolve_smf_context_from_tri)
    cov = df["count_meth"] + df["count_unmeth"]
    df["frac_meth"] = np.where(cov > 0, df["count_meth"] / cov.replace(0, 1), np.nan)
    df["pos0"] = df["pos"].astype(int) - 1  # 0-based for downstream code
    return df


def _resolve_smf_context_from_tri(tri: str) -> str:
    """Map a Bismark trinucleotide to an :class:`smf.context.Context` value (string)."""
    if not isinstance(tri, str) or len(tri) != 3:
        return Context.OTHER.value
    left, mid, right = tri.upper()
    if mid != "C":
        return Context.OTHER.value
    if left == "G" and right == "G":
        return Context.GCG.value
    if left == "G":
        return Context.GCH.value
    if right == "G":
        return Context.CCG.value if left == "C" else Context.HCG.value
    return Context.CHH.value


# ---------------------------------------------------------------------------
# Per-read methylation
# ---------------------------------------------------------------------------


@dataclass
class PerReadCall:
    """A single methylation call on a single read at a single cytosine."""

    read_id: str
    chrom: str
    pos: int           # 0-based reference position
    strand: str        # read strand: "+" or "-"
    context: Context
    methylated: bool   # True if methylated, False if unmethylated


def iter_xm_calls(
    bam_path: str | Path,
    region: tuple[str, int, int],
    fasta,
    *,
    min_mapq: int = 10,
) -> Iterator[PerReadCall]:
    """Yield :class:`PerReadCall`\\ s for every informative cytosine on every read
    overlapping ``region`` in the Bismark BAM at ``bam_path``.

    ``fasta`` should be a ``pyfaidx.Fasta`` (or any object with the same
    ``[chrom][start:end]`` slicing API returning a sequence-like with ``.seq``).

    This is the workhorse for building per-read matrices from real data; for the
    synthetic example we skip it and call :func:`smf.per_read.per_read_matrix` directly
    on a numpy array.
    """
    try:
        import pysam  # type: ignore[import-untyped]
    except ImportError as exc:  # pragma: no cover - exercised only with real data
        raise RuntimeError(
            "pysam is required to read BAMs; install it with `pip install pysam`."
        ) from exc

    chrom, start, end = region
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for read in bam.fetch(chrom, start, end):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.mapping_quality < min_mapq:
                continue
            try:
                xm = read.get_tag("XM")
            except KeyError:
                continue

            # Fetch one extra base on the left so that context_at() never hits its
            # pos==0 boundary guard when the detecting G is the very first base of
            # the aligned read (fixes GCH sites at the read's 5'-most position).
            fetch_start = max(0, read.reference_start - 1)
            ref_seq = str(fasta[chrom][fetch_start : read.reference_end + 1].seq).upper()

            # In Bismark non-directional paired-end BAM, whether XM characters appear
            # at C positions (use strand="+") or G positions (use strand="-") depends
            # on orientation, not simply on is_reverse:
            #   OT-R1 (FLAG≈99,  is_rev=F, is_r2=F): XM at C → "+"
            #   OT-R2 (FLAG≈147, is_rev=T, is_r2=T): XM at C → "+"
            #   OB-R1 (FLAG≈83,  is_rev=T, is_r2=F): XM at G → "-"
            #   OB-R2 (FLAG≈163, is_rev=F, is_r2=T): XM at G → "-"
            # Rule: strand="+" when is_reverse==is_read2, strand="-" otherwise.
            # Works for single-end too (is_read2=False, so is_reverse==False → "+").
            strand = "+" if (read.is_reverse == read.is_read2) else "-"

            ref_pos = read.reference_start
            qpos = 0
            for op, length in read.cigartuples or []:
                if op in (0, 7, 8):  # M / = / X -- consumes both
                    for k in range(length):
                        rp = ref_pos + k
                        offset = rp - fetch_start
                        if 0 <= offset < len(ref_seq):
                            ctx = context_at(ref_seq, offset, strand)
                            if ctx in (Context.GCH, Context.HCG, Context.CCG):
                                call_char = xm[qpos + k]
                                # Strand normalization: a duplex GpC site has its
                                # + strand C and its - strand C at DIFFERENT reference
                                # positions (one base apart). To make + and - strand
                                # calls of the same duplex site land in the same column
                                # of the per-read matrix, we always report the + strand
                                # C position regardless of which strand was sequenced.
                                if strand == "-":
                                    if ctx is Context.GCH:
                                        canonical_pos = rp + 1   # - strand C of GpC -> + strand C is one base downstream
                                    else:                         # HCG / CCG: both are CpG, - strand G at c+1 → + strand C at c
                                        canonical_pos = rp - 1   # - strand C of CpG -> + strand C is one base upstream
                                else:
                                    canonical_pos = rp
                                if call_char in "Zz" and ctx in (Context.HCG, Context.CCG):
                                    yield PerReadCall(
                                        read_id=read.query_name or "",
                                        chrom=chrom,
                                        pos=canonical_pos,
                                        strand=strand,
                                        context=Context.HCG,  # CCG is a CpG site too
                                        methylated=(call_char == "Z"),
                                    )
                                elif call_char in "HhXxUu" and ctx is Context.GCH:
                                    # In --CX mode Bismark uses the same letters for
                                    # any non-CpG context; SMF accessibility uses GCH.
                                    yield PerReadCall(
                                        read_id=read.query_name or "",
                                        chrom=chrom,
                                        pos=canonical_pos,
                                        strand=strand,
                                        context=Context.GCH,
                                        methylated=call_char.isupper(),
                                    )
                    ref_pos += length
                    qpos += length
                elif op == 1:  # I -- consumes query
                    qpos += length
                elif op == 2 or op == 3:  # D / N -- consumes ref
                    ref_pos += length
                elif op == 4:  # S -- consumes query
                    qpos += length
                # H (5), P (6) -- consume nothing


def read_per_read_tsv(
    path: str | Path,
    *,
    contexts: tuple[str, ...] = ("GCH",),
    chrom: str | None = None,
):
    """Read a per-read calls TSV (output of ``extract_per_read_methylation.py``)
    into a :class:`smf.per_read.PerReadMatrix`.

    The TSV is expected to have at least these columns (tab-separated, first row a
    header): ``read_id``, ``chrom``, ``pos``, ``strand``, ``context``, ``methylated``.

    Parameters
    ----------
    path
        Path to the TSV produced by the custom per-read caller.
    contexts
        Which contexts to keep -- usually ``("GCH",)`` for accessibility analysis,
        or ``("GCH", "HCG")`` if you want both layers.
    chrom
        If given, only keep calls on this chromosome (useful for multi-region TSVs).

    Returns
    -------
    smf.per_read.PerReadMatrix
        With one row per read and one column per informative position observed in
        the TSV. ``positions`` are 1-based, matching the TSV.
    """
    from .per_read import PerReadMatrix  # local import to avoid cycle

    df = pd.read_csv(path, sep="\t", dtype={"chrom": str, "strand": str, "context": str})
    if chrom is not None:
        df = df[df["chrom"] == chrom]
    df = df[df["context"].isin(contexts)].copy()
    if df.empty:
        raise ValueError(
            f"no rows match contexts={contexts}"
            + (f" on chrom={chrom!r}" if chrom is not None else "")
        )

    # Use a pivot to build the (n_reads x n_sites) matrix in one shot.
    df["pos"] = df["pos"].astype(int)
    df["methylated"] = df["methylated"].astype(int)
    pivot = df.pivot_table(
        index="read_id", columns="pos", values="methylated",
        aggfunc="mean",   # if a read covers a position twice (unlikely), average
    )
    sites = np.array(sorted(pivot.columns), dtype=int)
    pivot = pivot.reindex(columns=sites)
    matrix = pivot.to_numpy(dtype=float)
    read_ids = list(pivot.index)

    chrom_used = chrom if chrom is not None else (df["chrom"].iloc[0] if not df.empty else "")
    return PerReadMatrix(
        matrix=matrix,
        sites=sites,
        read_ids=read_ids,
        chrom=chrom_used,
        region_start=int(sites.min()) if len(sites) else 0,
        region_end=int(sites.max()) + 1 if len(sites) else 0,
    )


def calls_to_matrix(
    calls: Sequence[PerReadCall],
    sites: Sequence[int],
) -> tuple[np.ndarray, list[str]]:
    """Convert a flat list of :class:`PerReadCall` into a per-read methylation matrix.

    Parameters
    ----------
    calls
        Flat list of per-read calls (e.g. from :func:`iter_xm_calls`).
    sites
        Reference positions (0-based) defining the columns of the matrix, in order.

    Returns
    -------
    matrix
        ``float`` matrix of shape ``(n_reads, len(sites))`` with values in
        ``{0.0, 1.0, np.nan}`` (NaN = no call at that site for that read).
    read_ids
        The list of read IDs, indexing rows of ``matrix``.
    """
    site_index = {p: i for i, p in enumerate(sites)}

    by_read: dict[str, dict[int, bool]] = {}
    for c in calls:
        if c.pos in site_index:
            by_read.setdefault(c.read_id, {})[c.pos] = c.methylated

    read_ids = sorted(by_read)
    matrix = np.full((len(read_ids), len(sites)), np.nan, dtype=float)
    for r_idx, rid in enumerate(read_ids):
        for pos, meth in by_read[rid].items():
            matrix[r_idx, site_index[pos]] = 1.0 if meth else 0.0

    return matrix, read_ids
