#!/usr/bin/env python3
"""Write immutable-analysis provenance, final checksums, and completion marker."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import subprocess

import yaml


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def git_commit(root: Path) -> str:
    completed = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=root, check=False, capture_output=True, text=True
    )
    return completed.stdout.strip() if completed.returncode == 0 else "unavailable"


def configured_results_root(config: Path) -> Path:
    parsed = yaml.safe_load(config.read_text(encoding="utf-8")) or {}
    value = parsed.get("directories", {}).get("results", "")
    if not value:
        raise SystemExit(f"Missing directories.results in {config}")
    path = Path(value)
    if not path.is_absolute():
        path = config.resolve().parent / path
    return path.resolve()


def is_within(path: Path, root: Path) -> bool:
    return path == root or root in path.parents


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", required=True, type=Path)
    parser.add_argument("--result-dir", required=True, type=Path)
    parser.add_argument("--study", required=True)
    parser.add_argument("--analysis-id", required=True)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--metadata", action="append", type=Path, default=[])
    parser.add_argument("--container-lock", type=Path)
    parser.add_argument("--provenance", required=True, type=Path)
    parser.add_argument("--checksums", required=True, type=Path)
    parser.add_argument("--complete", required=True, type=Path)
    args = parser.parse_args()
    root = args.project_root.resolve()
    result_dir = args.result_dir.resolve()
    results_root = configured_results_root(args.config)
    if not (is_within(result_dir, root) or is_within(result_dir, results_root)):
        raise SystemExit(
            "Result directory must be inside project root or configured directories.results"
        )

    input_files = [args.config, *args.metadata]
    provenance = {
        "schema_version": "1.0",
        "study_id": args.study,
        "analysis_id": args.analysis_id,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "git_commit": git_commit(root),
        "platform": platform.platform(),
        "python": platform.python_version(),
        "inputs": {
            str(path): {"sha256": sha256(path), "bytes": path.stat().st_size}
            for path in input_files
        },
        "container_lock": (
            args.container_lock.read_text(encoding="utf-8") if args.container_lock and args.container_lock.is_file() else "missing"
        ),
    }
    args.provenance.parent.mkdir(parents=True, exist_ok=True)
    args.provenance.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    excluded = {args.checksums.resolve(), args.complete.resolve()}
    files = sorted(
        path for path in result_dir.rglob("*") if path.is_file() and path.resolve() not in excluded
    )
    with args.checksums.open("w", encoding="utf-8") as handle:
        for path in files:
            handle.write(f"{sha256(path)}  {path.relative_to(result_dir)}\n")
    args.complete.write_text(
        f"complete\nstudy={args.study}\nanalysis_id={args.analysis_id}\n", encoding="utf-8"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
