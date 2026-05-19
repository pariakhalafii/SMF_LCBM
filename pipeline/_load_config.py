"""Tiny YAML reader used by the shell pipeline.

Prints `KEY=VALUE` lines on stdout so they can be `eval`-ed by bash, e.g.:

    eval $(python pipeline/_load_config.py config/config.example.yaml)
    echo "$sample_id"

We only support flat scalar keys, which is all the pipeline needs.
"""
from __future__ import annotations

import sys
import shlex
from pathlib import Path


def _flat_yaml(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    for raw in path.read_text().splitlines():
        line = raw.split("#", 1)[0].rstrip()
        if not line.strip() or ":" not in line:
            continue
        key, _, val = line.partition(":")
        key = key.strip()
        val = val.strip().strip('"').strip("'")
        if not key or not val:
            continue
        out[key] = val
    return out


def main(argv: list[str]) -> int:
    if len(argv) != 2:
        print("usage: _load_config.py <config.yaml>", file=sys.stderr)
        return 2
    cfg = _flat_yaml(Path(argv[1]))
    for k, v in cfg.items():
        print(f"{k}={shlex.quote(v)}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
