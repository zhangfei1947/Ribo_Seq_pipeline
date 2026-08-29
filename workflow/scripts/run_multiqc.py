#!/usr/bin/env python3
"""Run MultiQC after exposing deterministically renamed FastQC archives."""

from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tempfile
import zipfile


def extract_fastqc_archives(search_paths: list[Path], destination: Path) -> int:
    extracted = 0
    for search_path in search_paths:
        if not search_path.exists():
            continue
        for archive in search_path.rglob("*.zip"):
            try:
                with zipfile.ZipFile(archive) as handle:
                    names = handle.namelist()
                    if not any(name.endswith("/fastqc_data.txt") for name in names):
                        continue
                    sample_dir = destination / f"archive_{extracted:04d}"
                    handle.extractall(sample_dir)
                    extracted += 1
            except zipfile.BadZipFile as error:
                raise SystemExit(f"Invalid ZIP archive: {archive}: {error}") from error
    return extracted


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--search-path", action="append", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--filename", default="multiqc_report.html")
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="multiqc-fastqc-", dir=args.outdir.parent) as tmp_name:
        extracted_dir = Path(tmp_name)
        extracted = extract_fastqc_archives(args.search_path, extracted_dir)
        command = ["multiqc", *(str(path) for path in args.search_path)]
        if extracted:
            command.append(str(extracted_dir))
        command.extend(["--force", "--outdir", str(args.outdir), "--filename", args.filename])
        subprocess.run(command, check=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
