#!/usr/bin/env python3
"""Retain primary, non-supplementary NH:i:1 alignments and report read fate."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pysam


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--metrics", required=True, type=Path)
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    counts = {"records": 0, "unmapped": 0, "secondary": 0, "supplementary": 0, "multimapping": 0, "unique": 0}
    with pysam.AlignmentFile(args.input, "rb") as source, pysam.AlignmentFile(
        args.output, "wb", template=source
    ) as destination:
        for read in source.fetch(until_eof=True):
            counts["records"] += 1
            if read.is_unmapped:
                counts["unmapped"] += 1
            elif read.is_secondary:
                counts["secondary"] += 1
            elif read.is_supplementary:
                counts["supplementary"] += 1
            elif not read.has_tag("NH") or read.get_tag("NH") != 1:
                counts["multimapping"] += 1
            else:
                destination.write(read)
                counts["unique"] += 1
    pysam.index(str(args.output))
    args.metrics.parent.mkdir(parents=True, exist_ok=True)
    args.metrics.write_text(json.dumps(counts, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

