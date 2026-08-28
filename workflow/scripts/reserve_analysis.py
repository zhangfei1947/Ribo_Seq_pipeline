#!/usr/bin/env python3
"""Prevent an analysis_id from silently acquiring different inputs or decisions."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path


def digest_inputs(paths: list[Path], study: str, analysis_id: str) -> str:
    digest = hashlib.sha256()
    digest.update(study.encode())
    digest.update(analysis_id.encode())
    for path in sorted(paths, key=str):
        digest.update(str(path).encode())
        digest.update(path.read_bytes())
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--result-dir", required=True, type=Path)
    parser.add_argument("--study", required=True)
    parser.add_argument("--analysis-id", required=True)
    parser.add_argument("--input", action="append", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    fingerprint = digest_inputs(args.input, args.study, args.analysis_id)
    result_dir = args.result_dir
    if result_dir.exists() and not args.output.exists():
        existing = [path for path in result_dir.iterdir() if path.name != args.output.name]
        if existing:
            raise SystemExit(
                f"Refusing to use non-empty result directory without analysis metadata: {result_dir}"
            )
    if args.output.exists():
        previous = json.loads(args.output.read_text(encoding="utf-8"))
        if previous.get("analysis_fingerprint") != fingerprint:
            raise SystemExit(
                f"analysis_id {args.analysis_id!r} already exists with different inputs; choose a new analysis_id"
            )
        return 0
    result_dir.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(
            {
                "schema_version": "1.0",
                "study_id": args.study,
                "analysis_id": args.analysis_id,
                "analysis_fingerprint": fingerprint,
                "reserved_utc": datetime.now(timezone.utc).isoformat(),
                "fingerprinted_inputs": [str(path) for path in sorted(args.input, key=str)],
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

