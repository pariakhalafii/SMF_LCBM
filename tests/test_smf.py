"""Tests for the smf package.

Run with:

    pytest -q

The tests do not require any real bioinformatics tools -- they exercise the package
on synthetic in-memory data.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

try:
    import pytest
except ImportError:  # pragma: no cover - allow tests to import without pytest
    class _MarkShim:
        def parametrize(self, argnames, argvalues):
            def decorator(func):
                # Match pytest's storage location so our run_tests.py can read it.
                marks = list(getattr(func, "pytestmark", []))
                marks.append(type("Mark", (), {"name": "parametrize", "args": (argnames, argvalues)})())
                func.pytestmark = marks
                return func
            return decorator

    class _PytestShim:
        mark = _MarkShim()

    pytest = _PytestShim()  # type: ignore[assignment]

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from scripts.make_synthetic_data import SyntheticConfig, make_synthetic  # noqa: E402
from smf import classify, context, per_read, stats  # noqa: E402


# ---------------------------------------------------------------------------
# context.py
# ---------------------------------------------------------------------------

@pytest.mark.parametrize(
    "tri,expected",
    [
        ("ACA", context.Context.CHH),
        ("ACG", context.Context.HCG),   # H = A here
        ("TCG", context.Context.HCG),
        ("CCG", context.Context.CCG),
        ("GCA", context.Context.GCH),
        ("GCT", context.Context.GCH),
        ("GCG", context.Context.GCG),   # ambiguous
        ("AAA", context.Context.OTHER),  # middle base not C
    ],
)
def test_classify_trinucleotide(tri, expected):
    assert context.classify_trinucleotide(tri) is expected


def test_revcomp_roundtrip():
    seq = "ACGTNGCATGCA"
    assert context.revcomp(context.revcomp(seq)) == seq


def test_iter_informative_sites_plus_and_minus():
    # ACGCATGCG -- positions 0-based:
    #  0 1 2 3 4 5 6 7 8
    #  A C G C A T G C G
    # Plus-strand C's are at positions 1, 3, 7.
    # Minus-strand C's (i.e. G on plus) are at 2, 6, 8.
    # Position 1 trinucleotide ACG -> HCG
    # Position 3 trinucleotide GCA -> GCH
    # Position 7 trinucleotide GCG -> GCG (excluded from default GCH-only iter)
    seq = "ACGCATGCG"
    sites = list(context.iter_informative_sites(seq, contexts=(context.Context.GCH,)))
    positions = sorted(s.pos for s in sites)
    # Expect at least position 3 on the plus strand.
    assert 3 in positions


# ---------------------------------------------------------------------------
# per_read.py
# ---------------------------------------------------------------------------

def test_per_read_matrix_dimensions_and_helpers():
    prm, _, _ = make_synthetic(SyntheticConfig(n_reads=50, seed=1))
    assert prm.matrix.shape == (prm.n_reads, prm.n_sites)
    assert len(prm.read_ids) == prm.n_reads
    cov = prm.coverage_per_read()
    assert cov.shape == (prm.n_reads,)
    assert np.all(cov <= prm.n_sites)


def test_filter_reads_removes_low_coverage():
    prm, _, _ = make_synthetic(SyntheticConfig(n_reads=200, seed=2))
    filt = prm.filter_reads(min_sites=10)
    assert filt.n_reads <= prm.n_reads
    assert np.all(filt.coverage_per_read() >= 10)


def test_sort_by_pattern_puts_protected_first():
    prm, labels, _ = make_synthetic(SyntheticConfig(n_reads=400, seed=3))
    prm = prm.filter_reads(min_sites=8)
    keep = np.isin(np.array([f"read_{i:05d}" for i in range(400)]), prm.read_ids)
    labels = labels[keep]

    sorted_prm = per_read.sort_by_pattern(prm, feature_center=1000, window=60)
    # Top quartile should be enriched for protected (TF-bound or nucleosome) molecules.
    top_q = max(1, sorted_prm.n_reads // 4)
    label_lookup = dict(zip(prm.read_ids, labels))
    top_labels = [label_lookup[r] for r in sorted_prm.read_ids[:top_q]]
    bottom_labels = [label_lookup[r] for r in sorted_prm.read_ids[-top_q:]]

    top_protected = sum(l in ("tf_bound", "nucleosome") for l in top_labels) / top_q
    bottom_protected = sum(l in ("tf_bound", "nucleosome") for l in bottom_labels) / top_q
    assert top_protected > bottom_protected


# ---------------------------------------------------------------------------
# classify.py
# ---------------------------------------------------------------------------

def test_classifier_recovers_majority_of_ground_truth():
    prm, true_labels, _ = make_synthetic(SyntheticConfig(n_reads=600, seed=4))
    prm_f = prm.filter_reads(min_sites=8)
    keep = np.isin(np.array([f"read_{i:05d}" for i in range(600)]), prm_f.read_ids)
    true_labels = true_labels[keep]

    result = classify.classify_footprints(prm_f, feature_center=1000)

    # For each ground-truth class, the classifier should call the correct state more
    # often than any other classified state.  We tolerate "unclassified" because some
    # synthetic reads have insufficient calls in the relevant window.
    for gt in ("accessible", "tf_bound", "nucleosome"):
        idx = true_labels == gt
        if idx.sum() < 10:
            continue
        preds = result.states[idx]
        classified = preds[preds != "unclassified"]
        if len(classified) == 0:
            continue
        unique, counts = np.unique(classified, return_counts=True)
        majority = unique[np.argmax(counts)]
        assert majority == gt, (
            f"Expected majority prediction {gt!r} for ground-truth {gt!r}, got {majority!r} "
            f"(distribution: {dict(zip(unique, counts))})"
        )


def test_classifier_handles_empty_window_gracefully():
    # Build a tiny matrix where the feature centre is outside every site.
    matrix = np.array([[1.0, 1.0, 1.0]])
    sites = np.array([0, 10, 20])
    prm = per_read.PerReadMatrix(matrix=matrix, sites=sites, read_ids=["r1"])
    result = classify.classify_footprints(prm, feature_center=10_000)
    assert list(result.states) == ["unclassified"]


# ---------------------------------------------------------------------------
# stats.py
# ---------------------------------------------------------------------------

def test_co_occupancy_perfect_correlation():
    # Construct a matrix where bin_a and bin_b are always equal -- log2(OR) should be > 0.
    sites = np.array([0, 1, 100, 101])
    rng = np.random.default_rng(0)
    n = 200
    a_vals = rng.integers(0, 2, size=n).astype(float)
    matrix = np.column_stack([a_vals, a_vals, a_vals, a_vals])
    prm = per_read.PerReadMatrix(matrix=matrix, sites=sites, read_ids=[f"r{i}" for i in range(n)])
    co = stats.co_occupancy(prm, bin_a=(0, 1), bin_b=(100, 101))
    assert co.n == n
    assert co.log2_or > 1.0


def test_average_accessibility_profile_returns_grid():
    prm, _, _ = make_synthetic(SyntheticConfig(n_reads=100, seed=5))
    grid, smoothed = stats.average_accessibility_profile(prm, smooth_bp=10)
    assert grid.shape == smoothed.shape
    assert np.all((smoothed >= 0) & (smoothed <= 1))
