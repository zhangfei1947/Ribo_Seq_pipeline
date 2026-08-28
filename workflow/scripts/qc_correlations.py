#!/usr/bin/env python3
"""Create pairwise replicate correlations from provisional gene count vectors."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def read_counts(path: Path) -> dict[str, float]:
    with path.open(encoding="utf-8") as handle:
        return {row["gene_id"]: math.log2(float(row["count"]) + 1.0) for row in csv.DictReader(handle, delimiter="\t")}


def pearson(left: list[float], right: list[float]) -> float | None:
    if len(left) < 2:
        return None
    left_mean, right_mean = sum(left) / len(left), sum(right) / len(right)
    numerator = sum((x - left_mean) * (y - right_mean) for x, y in zip(left, right))
    left_ss = sum((x - left_mean) ** 2 for x in left)
    right_ss = sum((y - right_mean) ** 2 for y in right)
    denominator = math.sqrt(left_ss * right_ss)
    return numerator / denominator if denominator else None


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--libraries", required=True, type=Path)
    parser.add_argument("--count", action="append", required=True, help="LIBRARY=PATH")
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    with args.libraries.open(encoding="utf-8") as handle:
        metadata = {row["library_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
    counts = {item.split("=", 1)[0]: read_counts(Path(item.split("=", 1)[1])) for item in args.count}
    rows: list[dict[str, object]] = []
    library_ids = sorted(counts)
    for index, left_id in enumerate(library_ids):
        for right_id in library_ids[index + 1 :]:
            left_meta, right_meta = metadata[left_id], metadata[right_id]
            if (left_meta["study_id"], left_meta["assay"], left_meta["condition"]) != (
                right_meta["study_id"], right_meta["assay"], right_meta["condition"]
            ):
                continue
            genes = sorted(set(counts[left_id]) & set(counts[right_id]))
            value = pearson([counts[left_id][gene] for gene in genes], [counts[right_id][gene] for gene in genes])
            rows.append(
                {
                    "study_id": left_meta["study_id"],
                    "assay": left_meta["assay"],
                    "condition": left_meta["condition"],
                    "library_1": left_id,
                    "library_2": right_id,
                    "pearson_log2_count_plus_1": "" if value is None else f"{value:.6f}",
                    "n_genes": len(genes),
                }
            )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fields = ["study_id", "assay", "condition", "library_1", "library_2", "pearson_log2_count_plus_1", "n_genes"]
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

