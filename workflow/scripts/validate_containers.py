#!/usr/bin/env python3
"""Verify production SIF files against the committed container lock."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--lock", required=True, type=Path)
    parser.add_argument("--image", action="append", required=True, help="ROLE=PATH")
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    expected = {item.split("=", 1)[0]: Path(item.split("=", 1)[1]) for item in args.image}
    with args.lock.open(encoding="utf-8") as handle:
        rows = {row["role"]: row for row in csv.DictReader(handle, delimiter="\t")}
    report: dict[str, object] = {"status": "valid", "images": {}}
    for role, path in expected.items():
        if role not in rows:
            raise SystemExit(f"Container lock is missing role {role}")
        if not path.is_file():
            raise SystemExit(f"Missing SIF for role {role}: {path}")
        observed = sha256(path)
        locked = rows[role]["sha256"].lower()
        if len(locked) != 64 or observed != locked:
            raise SystemExit(f"SIF checksum mismatch for {role}: expected {locked}, observed {observed}")
        report["images"][role] = {"path": str(path), "sha256": observed, "software": rows[role]["software"]}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

