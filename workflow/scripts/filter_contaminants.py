#!/usr/bin/env python3
"""Remove Ribo reads with a credible end-to-end contaminant alignment."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--index", required=True)
    parser.add_argument("--threads", default=1, type=int)
    parser.add_argument("--metrics", required=True, type=Path)
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.metrics.parent.mkdir(parents=True, exist_ok=True)
    command = [
        "bowtie2",
        "--very-sensitive",
        "--end-to-end",
        "--threads",
        str(args.threads),
        "-x",
        args.index,
        "-U",
        str(args.input),
        "--un-gz",
        str(args.output),
        "-S",
        "/dev/null",
    ]
    completed = subprocess.run(command, check=True, capture_output=True, text=True)
    args.metrics.write_text(
        json.dumps(
            {"status": "complete", "policy": "any_reported_end_to_end_hit_removed", "bowtie2_stderr": completed.stderr},
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

