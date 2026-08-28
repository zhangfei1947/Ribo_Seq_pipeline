#!/usr/bin/env python3
"""Create a review template and a fingerprint bound to all QC inputs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path


FIELDS = [
    "library_id",
    "include_assay_analysis",
    "include_te",
    "selected_lengths",
    "offsets",
    "offset_source",
    "adapter_review",
    "confirmed_adapter_3p",
    "qc_status",
    "reason",
    "reviewer",
    "review_date",
    "decision_status",
    "qc_fingerprint",
]


def fingerprint(paths: list[Path]) -> str:
    digest = hashlib.sha256()
    for path in sorted(paths, key=str):
        digest.update(str(path).encode())
        digest.update(path.read_bytes())
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--libraries", required=True, type=Path)
    parser.add_argument("--qc-files", nargs="*", type=Path, default=[])
    parser.add_argument("--offset-file", action="append", default=[], help="LIBRARY=PATH")
    parser.add_argument("--pooled-offset-file", action="append", default=[], help="LIBRARY=PATH")
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--fingerprint", required=True, type=Path)
    args = parser.parse_args()
    digest = fingerprint([args.libraries, *args.qc_files])
    offset_paths = {item.split("=", 1)[0]: Path(item.split("=", 1)[1]) for item in args.offset_file}
    pooled_paths = {item.split("=", 1)[0]: Path(item.split("=", 1)[1]) for item in args.pooled_offset_file}
    with args.libraries.open(encoding="utf-8") as handle:
        libraries = list(csv.DictReader(handle, delimiter="\t"))
    rows: list[dict[str, str]] = []
    for library in libraries:
        selected_lengths = ""
        offsets = ""
        offset_source = ""
        if library["assay"] == "ribo" and library["library_id"] in offset_paths:
            with offset_paths[library["library_id"]].open(encoding="utf-8") as handle:
                recommendations = [
                    row for row in csv.DictReader(handle, delimiter="\t") if row["status"] == "recommended"
                ]
            if not recommendations and library["library_id"] in pooled_paths:
                with pooled_paths[library["library_id"]].open(encoding="utf-8") as handle:
                    recommendations = [
                        row for row in csv.DictReader(handle, delimiter="\t") if row["status"] == "recommended"
                    ]
                offset_source = "pooled"
            selected_lengths = ",".join(row["read_length"] for row in recommendations)
            offsets = ",".join(
                f"{row['read_length']}:{row['recommended_offset']}" for row in recommendations
            )
            if not offset_source:
                offset_source = "sample"
        adapter = library.get("adapter_3p", "")
        adapter_review = "REVIEW_REQUIRED" if adapter == "infer" else "none" if adapter in {"", "none"} else "explicit"
        confirmed_adapter = "" if adapter == "infer" else adapter
        rows.append(
            {
                "library_id": library["library_id"],
                "include_assay_analysis": "REVIEW_REQUIRED",
                "include_te": "REVIEW_REQUIRED",
                "selected_lengths": selected_lengths,
                "offsets": offsets,
                "offset_source": offset_source,
                "adapter_review": adapter_review,
                "confirmed_adapter_3p": confirmed_adapter,
                "qc_status": "REVIEW_REQUIRED",
                "reason": "",
                "reviewer": "",
                "review_date": "",
                "decision_status": "pending",
                "qc_fingerprint": digest,
            }
        )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    args.fingerprint.parent.mkdir(parents=True, exist_ok=True)
    args.fingerprint.write_text(digest + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
