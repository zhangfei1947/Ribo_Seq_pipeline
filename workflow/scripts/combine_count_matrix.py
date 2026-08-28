#!/usr/bin/env python3
"""Join per-library integer counts with stable gene_id/gene_symbol annotation."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_counts(path: Path) -> dict[str, int]:
    with path.open(encoding="utf-8") as handle:
        return {row["gene_id"]: int(row["count"]) for row in csv.DictReader(handle, delimiter="\t")}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--annotation", required=True, type=Path)
    parser.add_argument("--count", action="append", required=True, help="LIBRARY=PATH")
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    entries = [item.split("=", 1) for item in args.count]
    if len({library for library, _ in entries}) != len(entries):
        raise SystemExit("Duplicate library in --count inputs")
    matrices = {library: read_counts(Path(path)) for library, path in entries}
    with args.annotation.open(encoding="utf-8") as handle:
        annotation = list(csv.DictReader(handle, delimiter="\t"))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fields = ["gene_id", "gene_symbol", *[library for library, _ in entries]]
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for gene in annotation:
            gene_id = gene["gene_id"]
            row: dict[str, object] = {"gene_id": gene_id, "gene_symbol": gene["gene_symbol"]}
            row.update({library: counts.get(gene_id, 0) for library, counts in matrices.items()})
            writer.writerow(row)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

