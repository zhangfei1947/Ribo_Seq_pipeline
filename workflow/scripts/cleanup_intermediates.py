#!/usr/bin/env python3
"""Safely delete only exact, manifest-listed large intermediates; dry-run by default."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


ALLOWED_SUBDIRECTORIES = {"preprocessed", "libraries", "depleted", "alignment"}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", default=".", type=Path)
    parser.add_argument("--work-dir", default="work", type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--analysis-id", required=True)
    parser.add_argument("--execute", action="store_true")
    parser.add_argument("--confirm-analysis-id", default="")
    args = parser.parse_args()
    root = args.project_root.resolve()
    work = (root / args.work_dir).resolve() if not args.work_dir.is_absolute() else args.work_dir.resolve()
    if root not in work.parents:
        raise SystemExit("Work directory must be inside project root")
    if args.execute and args.confirm_analysis_id != args.analysis_id:
        raise SystemExit("Execution requires --confirm-analysis-id equal to --analysis-id")

    with args.manifest.open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    total = 0
    targets: list[Path] = []
    for row in rows:
        if row.get("analysis_id") != args.analysis_id:
            raise SystemExit("Manifest analysis_id mismatch")
        relative = Path(row["relative_path"])
        if relative.is_absolute() or ".." in relative.parts:
            raise SystemExit(f"Unsafe manifest path: {relative}")
        target = (root / relative).resolve()
        if work not in target.parents:
            raise SystemExit(f"Target is outside work directory: {target}")
        under_work = target.relative_to(work)
        if not under_work.parts or under_work.parts[0] not in ALLOWED_SUBDIRECTORIES:
            raise SystemExit(f"Disallowed cleanup category: {target}")
        if target.is_dir():
            raise SystemExit(f"Manifest may contain files only: {target}")
        if target.exists() or target.is_symlink():
            total += target.lstat().st_size
            targets.append(target)

    action = "DELETE" if args.execute else "WOULD_DELETE"
    for target in targets:
        print(f"{action}\t{target}")
        if args.execute:
            target.unlink()
    print(f"files={len(targets)} bytes={total} mode={'execute' if args.execute else 'dry-run'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

