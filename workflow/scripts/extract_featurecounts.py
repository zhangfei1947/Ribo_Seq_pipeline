#!/usr/bin/env python3
"""Convert a featureCounts table to a stable two-column gene count TSV."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.input.open(encoding="utf-8") as source, args.output.open(
        "w", newline="", encoding="utf-8"
    ) as destination:
        rows = (line for line in source if not line.startswith("#"))
        reader = csv.DictReader(rows, delimiter="\t")
        if not reader.fieldnames or "Geneid" not in reader.fieldnames:
            raise SystemExit(f"Invalid featureCounts file: {args.input}")
        count_column = reader.fieldnames[-1]
        writer = csv.writer(destination, delimiter="\t", lineterminator="\n")
        writer.writerow(["gene_id", "count"])
        for row in reader:
            writer.writerow([row["Geneid"], int(float(row[count_column]))])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

