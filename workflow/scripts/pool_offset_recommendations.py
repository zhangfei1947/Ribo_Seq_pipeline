#!/usr/bin/env python3
"""Pool candidate frame evidence within a study/protocol group."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
from pathlib import Path


FIELDS = [
    "read_length", "recommended_offset", "n_representative_cds", "frame0", "frame1", "frame2",
    "frame0_fraction", "second_best_margin", "status",
]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate", action="append", required=True, type=Path)
    parser.add_argument("--min-reads-per-length", required=True, type=int)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    pooled: dict[tuple[int, int], Counter[int]] = defaultdict(Counter)
    for path in args.candidate:
        with path.open(encoding="utf-8") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                key = (int(row["read_length"]), int(row["offset"]))
                for frame in range(3):
                    pooled[key][frame] += int(row[f"frame{frame}"])
    by_length: dict[int, list[tuple[float, int, int, Counter[int]]]] = defaultdict(list)
    for (length, offset), frames in pooled.items():
        total = sum(frames.values())
        score = frames[0] / total if total else 0.0
        by_length[length].append((score, total, -offset, frames))
    rows: list[dict[str, object]] = []
    for length in sorted(by_length):
        ranked = sorted(by_length[length], reverse=True)
        score, total, negative_offset, frames = ranked[0]
        second = ranked[1][0] if len(ranked) > 1 else 0.0
        rows.append(
            {
                "read_length": length,
                "recommended_offset": -negative_offset,
                "n_representative_cds": total,
                "frame0": frames[0],
                "frame1": frames[1],
                "frame2": frames[2],
                "frame0_fraction": f"{score:.6f}",
                "second_best_margin": f"{score - second:.6f}",
                "status": "recommended" if total >= args.min_reads_per_length else "insufficient_reads",
            }
        )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
