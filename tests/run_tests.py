"""Tiny stdlib-only test runner for environments without pytest installed.

Discovers every ``test_*`` function in ``tests/test_smf.py`` (including parametrised
ones) and runs them, printing pass/fail per test.

For real development use ``pytest -q`` instead -- it gives much better tracebacks.
"""
from __future__ import annotations

import importlib
import inspect
import sys
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))


def _expand_parametrize(func):
    """Return a list of (label, callable) pairs, expanding pytest.mark.parametrize."""
    params = getattr(func, "pytestmark", [])
    cases = []
    for mark in params:
        if getattr(mark, "name", None) == "parametrize":
            arg_names = mark.args[0]
            arg_values = mark.args[1]
            names = [n.strip() for n in arg_names.split(",")]
            for vals in arg_values:
                if not isinstance(vals, tuple):
                    vals = (vals,)
                kwargs = dict(zip(names, vals))
                label = ",".join(f"{k}={v!r}" for k, v in kwargs.items())
                cases.append((label, lambda kw=kwargs: func(**kw)))
            return cases
    if cases:
        return cases
    return [("", func)]


def main() -> int:
    mod = importlib.import_module("tests.test_smf")
    funcs = [
        (name, obj)
        for name, obj in inspect.getmembers(mod, inspect.isfunction)
        if name.startswith("test_")
    ]
    funcs.sort(key=lambda kv: kv[0])

    passed = failed = 0
    for name, func in funcs:
        for label, runner in _expand_parametrize(func):
            full = f"{name}[{label}]" if label else name
            try:
                runner()
                print(f"  PASS  {full}")
                passed += 1
            except Exception:
                print(f"  FAIL  {full}")
                traceback.print_exc()
                failed += 1
    print(f"\n{passed} passed, {failed} failed")
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
