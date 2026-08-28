#!/usr/bin/env python3
"""Enumerate exact disposable files; never place counts or final results in the manifest."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


ALLOWED_SUBDIRECTORIES = ("preprocessed", "libraries", "depleted", "alignment")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", required=True, type=Path)
    parser.add_argument("--work-dir", required=True, type=Path)
    parser.add_argument("--analysis-id", required=True)
    parser.add_argument("--completion", action="append", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    root = args.project_root.resolve()
    work = (root / args.work_dir).resolve() if not args.work_dir.is_absolute() else args.work_dir.resolve()
    for marker in args.completion:
        if not marker.is_file() or not marker.read_text(encoding="utf-8").startswith("complete\n"):
            raise SystemExit(f"Analysis is not finalized: {marker}")

    rows: list[dict[str, object]] = []
    for subdirectory in ALLOWED_SUBDIRECTORIES:
        directory = work / subdirectory
        if not directory.is_dir():
            continue
        for path in sorted(directory.rglob("*")):
            if path.is_file() or path.is_symlink():
                rows.append(
                    {
                        "relative_path": str(path.resolve().relative_to(root)),
                        "bytes": path.stat().st_size,
                        "category": subdirectory,
                        "analysis_id": args.analysis_id,
                    }
                )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["relative_path", "bytes", "category", "analysis_id"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

