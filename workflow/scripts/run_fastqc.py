#!/usr/bin/env python3
"""Run FastQC while giving Snakemake deterministic output names."""

from __future__ import annotations

import argparse
from pathlib import Path
import shutil
import subprocess
import tempfile


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--html", required=True, type=Path)
    parser.add_argument("--zip", required=True, type=Path)
    parser.add_argument("--threads", default=1, type=int)
    args = parser.parse_args()
    args.html.parent.mkdir(parents=True, exist_ok=True)
    args.zip.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="fastqc-", dir=args.html.parent) as tmp_name:
        subprocess.run(
            ["fastqc", "--threads", str(args.threads), "--outdir", tmp_name, str(args.input)],
            check=True,
        )
        html_files = list(Path(tmp_name).glob("*_fastqc.html"))
        zip_files = list(Path(tmp_name).glob("*_fastqc.zip"))
        if len(html_files) != 1 or len(zip_files) != 1:
            raise SystemExit("FastQC did not create exactly one HTML and ZIP output")
        shutil.move(html_files[0], args.html)
        shutil.move(zip_files[0], args.zip)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

