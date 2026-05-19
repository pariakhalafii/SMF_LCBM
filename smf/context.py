"""Sequence-context utilities for NOMe-seq / SMF.

In dual-enzyme bisulfite SMF, a cytosine's biological meaning depends entirely on its
*trinucleotide context* on the **plus** strand:

* **GCH** (G-C-H, where H = A/C/T): a **GpC** site that is *not* also a CpG.
  Methylation here is from the M.CviPI treatment and reports **accessibility**.
* **HCG** (H-C-G): a **CpG** site that is *not* also a GpC. Methylation here is the
  cell's **endogenous** methylation.
* **GCG**: ambiguous -- the C is part of *both* a GpC and a CpG context. These sites
  cannot be unambiguously assigned to either signal and are excluded from analysis.
* **CCG / CHH / CHG**: not informative for SMF; ignored.

For positions on the minus strand the context is determined from the reverse-complement,
i.e. the C is the central base of the trinucleotide ``rev_comp(seq[pos-1:pos+2])``.

This module contains pure-Python helpers; reading FASTA is delegated to ``pyfaidx`` if
available, but everything works on bare strings too so it is trivially testable.
"""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Iterable, Iterator


class Context(str, Enum):
    """Resolved SMF context for a cytosine."""

    GCH = "GCH"   # informative GpC -- accessibility readout
    HCG = "HCG"   # informative CpG -- endogenous methylation
    GCG = "GCG"   # ambiguous -- exclude
    CCG = "CCG"   # other -- ignore
    CHH = "CHH"   # other -- ignore
    CHG = "CHG"   # other -- ignore
    OTHER = "OTHER"


_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq: str) -> str:
    """Reverse complement of a DNA string."""
    return seq.translate(_COMPLEMENT)[::-1]


def classify_trinucleotide(tri: str) -> Context:
    """Classify a 3-letter trinucleotide centred on a C.

    Parameters
    ----------
    tri : str
        Three uppercase DNA bases. The middle base must be ``C``; otherwise
        :data:`Context.OTHER` is returned.
    """
    if len(tri) != 3 or tri[1] != "C":
        return Context.OTHER

    left, _, right = tri
    is_gpc = left == "G"
    is_cpg = right == "G"

    if is_gpc and is_cpg:
        return Context.GCG
    if is_gpc:
        return Context.GCH
    if is_cpg:
        # Distinguish HCG (informative CpG) from CCG (treated separately by Bismark).
        if left == "C":
            return Context.CCG
        return Context.HCG
    if right == "G":
        return Context.CHG  # unreachable (would be CpG above) but kept for completeness
    return Context.CHH


def context_at(seq: str, pos: int, strand: str = "+") -> Context:
    """Return the SMF context of position ``pos`` (0-based) within ``seq``.

    For strand ``"-"`` the trinucleotide is reverse-complemented before classification,
    matching how Bismark reports CX contexts.
    """
    if pos <= 0 or pos >= len(seq) - 1:
        return Context.OTHER
    tri = seq[pos - 1 : pos + 2].upper()
    if strand == "-":
        tri = revcomp(tri)
    return classify_trinucleotide(tri)


@dataclass(frozen=True)
class CytosineSite:
    """A single informative cytosine in some reference window."""

    pos: int        # 0-based position within the window
    context: Context
    strand: str     # "+" or "-"


def iter_informative_sites(
    seq: str,
    *,
    contexts: Iterable[Context] = (Context.GCH,),
    include_minus: bool = True,
) -> Iterator[CytosineSite]:
    """Yield :class:`CytosineSite` for every informative cytosine in ``seq``.

    By default only ``GCH`` (accessibility) sites are yielded; pass
    ``contexts=(Context.GCH, Context.HCG)`` to include endogenous CpGs as well.
    """
    wanted = set(contexts)
    upper = seq.upper()
    for i, base in enumerate(upper):
        if base == "C":
            ctx = context_at(upper, i, "+")
            if ctx in wanted:
                yield CytosineSite(pos=i, context=ctx, strand="+")
        elif include_minus and base == "G":
            ctx = context_at(upper, i, "-")
            if ctx in wanted:
                yield CytosineSite(pos=i, context=ctx, strand="-")
