#!/usr/bin/env python3
"""Run STAR with frozen assay-specific alignment policy."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil
import subprocess
import tempfile


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--read1", required=True, type=Path)
    parser.add_argument("--read2", default="", type=Path)
    parser.add_argument("--layout", required=True, choices=("SE", "PE"))
    parser.add_argument("--assay", required=True, choices=("ribo", "rna"))
    parser.add_argument("--index", required=True, type=Path)
    parser.add_argument("--output-bam", required=True, type=Path)
    parser.add_argument("--log-final", required=True, type=Path)
    parser.add_argument("--metrics", required=True, type=Path)
    parser.add_argument("--threads", default=1, type=int)
    parser.add_argument("--mismatches", default=2, type=int)
    parser.add_argument("--multimap-max", default=20, type=int)
    args = parser.parse_args()

    for path in (args.output_bam, args.log_final, args.metrics):
        path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="star-", dir=args.output_bam.parent) as tmp_name:
        prefix = str(Path(tmp_name) / "star.")
        command = [
            "STAR",
            "--runThreadN",
            str(args.threads),
            "--genomeDir",
            str(args.index),
            "--readFilesIn",
            str(args.read1),
        ]
        if args.layout == "PE":
            command.append(str(args.read2))
        command.extend(
            [
                "--readFilesCommand",
                "zcat",
                "--outFileNamePrefix",
                prefix,
                "--outSAMtype",
                "BAM",
                "SortedByCoordinate",
                "--outSAMattributes",
                "NH",
                "HI",
                "AS",
                "nM",
                "--outFilterMultimapNmax",
                str(args.multimap_max),
            ]
        )
        if args.assay == "ribo":
            command.extend(
                [
                    "--alignEndsType",
                    "EndToEnd",
                    "--outFilterMismatchNmax",
                    str(args.mismatches),
                    "--outFilterMismatchNoverLmax",
                    "1",
                ]
            )
        subprocess.run(command, check=True)
        shutil.move(prefix + "Aligned.sortedByCoord.out.bam", args.output_bam)
        shutil.copy2(prefix + "Log.final.out", args.log_final)
    subprocess.run(["samtools", "index", "-@", str(args.threads), str(args.output_bam)], check=True)
    args.metrics.write_text(
        json.dumps(
            {
                "status": "complete",
                "assay": args.assay,
                "layout": args.layout,
                "unique_policy_downstream": "primary_non_supplementary_NH_eq_1",
                "ribo_end_to_end": args.assay == "ribo",
                "mismatch_max": args.mismatches if args.assay == "ribo" else None,
                "multimap_report_max": args.multimap_max,
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

