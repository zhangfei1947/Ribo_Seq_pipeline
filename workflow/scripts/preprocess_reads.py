#!/usr/bin/env python3
"""Parameter-driven UMI extraction, adapter trimming, and 3' quality trimming."""

from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path
import shutil
import subprocess
import tempfile


def run(command: list[str]) -> None:
    subprocess.run(command, check=True)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fastq1", required=True, type=Path)
    parser.add_argument("--fastq2", default="", type=Path)
    parser.add_argument("--output1", required=True, type=Path)
    parser.add_argument("--output2", required=True, type=Path)
    parser.add_argument("--layout", required=True, choices=("SE", "PE"))
    parser.add_argument("--adapter", required=True)
    parser.add_argument("--umi-pattern", default="none")
    parser.add_argument("--minimum-length", required=True, type=int)
    parser.add_argument("--maximum-length", type=int)
    parser.add_argument("--quality-3p", default=0, type=int)
    parser.add_argument("--threads", default=1, type=int)
    parser.add_argument("--report", required=True, type=Path)
    args = parser.parse_args()

    if args.layout == "PE" and not str(args.fastq2):
        raise SystemExit("PE preprocessing requires --fastq2")
    args.output1.parent.mkdir(parents=True, exist_ok=True)
    args.output2.parent.mkdir(parents=True, exist_ok=True)
    args.report.parent.mkdir(parents=True, exist_ok=True)

    report: dict[str, object] = {
        "layout": args.layout,
        "adapter_requested": args.adapter,
        "umi_pattern": args.umi_pattern,
        "minimum_length": args.minimum_length,
        "maximum_length": args.maximum_length,
        "quality_3p": args.quality_3p,
    }
    with tempfile.TemporaryDirectory(prefix="preprocess-", dir=args.output1.parent) as tmp_name:
        tmp = Path(tmp_name)
        input1 = args.fastq1
        input2 = args.fastq2
        if args.umi_pattern not in {"", "none", "NA", "."}:
            umi1 = tmp / "umi_R1.fastq.gz"
            umi2 = tmp / "umi_R2.fastq.gz"
            command = [
                "umi_tools",
                "extract",
                "--bc-pattern",
                args.umi_pattern,
                "--stdin",
                str(input1),
                "--stdout",
                str(umi1),
            ]
            if args.layout == "PE":
                command.extend(["--read2-in", str(input2), "--read2-out", str(umi2)])
            run(command)
            input1, input2 = umi1, umi2

        if args.adapter == "infer":
            # Inference is recorded and must later be accepted in the reviewed QC file.
            command = [
                "fastp",
                "--in1",
                str(input1),
                "--out1",
                str(args.output1),
                "--thread",
                str(args.threads),
                "--length_required",
                str(args.minimum_length),
                "--json",
                str(tmp / "fastp.json"),
                "--html",
                str(tmp / "fastp.html"),
            ]
            if args.maximum_length:
                command.extend(["--length_limit", str(args.maximum_length)])
            if args.quality_3p:
                command.extend(["--cut_right", "--cut_right_mean_quality", str(args.quality_3p)])
            if args.layout == "PE":
                command.extend(
                    ["--in2", str(input2), "--out2", str(args.output2), "--detect_adapter_for_pe"]
                )
            run(command)
            fastp_report = json.loads((tmp / "fastp.json").read_text(encoding="utf-8"))
            report["adapter_method"] = "fastp_inference_requires_review"
            report["fastp_adapter_cutting"] = fastp_report.get("adapter_cutting", {})
        else:
            command = ["cutadapt", "-j", str(args.threads), "-m", str(args.minimum_length)]
            if args.maximum_length:
                command.extend(["-M", str(args.maximum_length)])
            if args.quality_3p:
                command.extend(["-q", f"0,{args.quality_3p}"])
            if args.adapter not in {"none", "NA", ".", ""}:
                command.extend(["-a", args.adapter])
                if args.layout == "PE":
                    command.extend(["-A", args.adapter])
                report["adapter_method"] = "explicit"
            else:
                report["adapter_method"] = "none"
            command.extend(["-o", str(args.output1)])
            if args.layout == "PE":
                command.extend(["-p", str(args.output2), str(input1), str(input2)])
            else:
                command.append(str(input1))
            run(command)

    if args.layout == "SE":
        with gzip.open(args.output2, "wb"):
            pass
    report["status"] = "complete"
    args.report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

