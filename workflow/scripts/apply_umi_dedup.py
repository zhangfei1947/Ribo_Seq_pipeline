#!/usr/bin/env python3
"""Apply UMI-aware deduplication only when the metadata declares a UMI pattern."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil
import subprocess


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--input-index", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--output-index", required=True, type=Path)
    parser.add_argument("--umi-pattern", default="none")
    parser.add_argument("--layout", required=True, choices=("SE", "PE"))
    parser.add_argument("--metrics", required=True, type=Path)
    parser.add_argument("--log", required=True, type=Path)
    args = parser.parse_args()
    for path in (args.output, args.output_index, args.metrics, args.log):
        path.parent.mkdir(parents=True, exist_ok=True)
    enabled = args.umi_pattern not in {"", "none", "NA", "."}
    if enabled:
        command = [
            "umi_tools",
            "dedup",
            "--stdin",
            str(args.input),
            "--stdout",
            str(args.output),
            "--log",
            str(args.log),
        ]
        if args.layout == "PE":
            command.append("--paired")
        subprocess.run(command, check=True)
        subprocess.run(["samtools", "index", str(args.output)], check=True)
    else:
        shutil.copy2(args.input, args.output)
        shutil.copy2(args.input_index, args.output_index)
        args.log.write_text("UMI deduplication disabled: no explicit UMI pattern\n", encoding="utf-8")
    args.metrics.write_text(
        json.dumps(
            {
                "umi_deduplication": "enabled" if enabled else "disabled",
                "coordinate_only_deduplication": False,
                "umi_pattern": args.umi_pattern,
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

