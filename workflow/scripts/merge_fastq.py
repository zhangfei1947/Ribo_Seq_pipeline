#!/usr/bin/env python3
"""Concatenate FASTQ gzip members without recompressing them."""

from __future__ import annotations

import argparse
from pathlib import Path
import shutil


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputs", nargs="+", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("wb") as destination:
        for source_path in args.inputs:
            with source_path.open("rb") as source:
                shutil.copyfileobj(source, destination, length=1024 * 1024)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

