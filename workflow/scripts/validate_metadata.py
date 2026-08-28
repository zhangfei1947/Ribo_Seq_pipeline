#!/usr/bin/env python3
"""Strictly validate sample, contrast, and optional reviewed-QC metadata."""

from __future__ import annotations

import argparse
from collections import defaultdict
import csv
from datetime import date
import hashlib
import json
from pathlib import Path
import re
import sys
from typing import Iterable


SAMPLE_REQUIRED = {
    "study_id",
    "run_id",
    "library_id",
    "bio_sample_id",
    "assay",
    "condition",
    "fastq1",
    "layout",
}
CONTRAST_REQUIRED = {
    "study_id",
    "contrast_id",
    "test_condition",
    "reference_condition",
    "run_ribo_de",
    "run_rna_de",
    "run_dte",
}
DECISION_REQUIRED = {
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
}
YES_NO = {"yes", "no"}
SAFE_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")


class ValidationError(ValueError):
    pass


def read_tsv(path: Path, required: set[str]) -> tuple[list[str], list[dict[str, str]]]:
    if not path.is_file():
        raise ValidationError(f"Missing TSV: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        missing = sorted(required - set(fields))
        if missing:
            raise ValidationError(f"{path}: missing columns: {', '.join(missing)}")
        rows = []
        for line_number, row in enumerate(reader, start=2):
            cleaned = {key: (value or "").strip() for key, value in row.items()}
            if not any(cleaned.values()):
                continue
            cleaned["_line"] = str(line_number)
            rows.append(cleaned)
    if not rows:
        raise ValidationError(f"{path}: contains no data rows")
    return fields, rows


def require_value(row: dict[str, str], field: str, table: Path) -> str:
    value = row.get(field, "").strip()
    if not value:
        raise ValidationError(f"{table}:{row['_line']}: {field} is empty")
    return value


def consistent(rows: Iterable[dict[str, str]], field: str, label: str) -> str:
    values = {row.get(field, "") for row in rows}
    if len(values) != 1:
        raise ValidationError(f"{label} has inconsistent {field}: {sorted(values)}")
    return next(iter(values))


def validate_samples(path: Path, check_fastq: bool) -> tuple[list[dict[str, str]], list[dict[str, str]], list[str]]:
    fields, rows = read_tsv(path, SAMPLE_REQUIRED)
    errors: list[str] = []
    seen_runs: set[str] = set()
    by_library: dict[str, list[dict[str, str]]] = defaultdict(list)

    for row in rows:
        for field in SAMPLE_REQUIRED:
            require_value(row, field, path)
        for field in ("study_id", "run_id", "library_id", "bio_sample_id"):
            if not SAFE_ID.fullmatch(row[field]):
                raise ValidationError(
                    f"{path}:{row['_line']}: {field} must match {SAFE_ID.pattern}"
                )
        run_id = row["run_id"]
        if run_id in seen_runs:
            raise ValidationError(f"{path}:{row['_line']}: duplicate run_id {run_id}")
        seen_runs.add(run_id)
        if row["assay"] not in {"ribo", "rna"}:
            raise ValidationError(f"{path}:{row['_line']}: assay must be ribo or rna")
        if row["layout"] not in {"SE", "PE"}:
            raise ValidationError(f"{path}:{row['_line']}: layout must be SE or PE")
        if row["assay"] == "ribo" and row["layout"] != "SE":
            raise ValidationError(f"{path}:{row['_line']}: v1 supports single-end Ribo-seq only")
        if row.get("strandedness", "") not in {"forward", "reverse", "unstranded"}:
            raise ValidationError(
                f"{path}:{row['_line']}: strandedness must be forward, reverse, or unstranded"
            )
        if row["assay"] == "ribo" and row.get("strandedness") == "unstranded":
            raise ValidationError(f"{path}:{row['_line']}: v1 Ribo-seq counting requires stranded data")
        if row["assay"] == "ribo" and not row.get("ribo_protocol_group", ""):
            raise ValidationError(f"{path}:{row['_line']}: Ribo run requires ribo_protocol_group")
        if row["layout"] == "PE" and not row.get("fastq2", ""):
            raise ValidationError(f"{path}:{row['_line']}: PE run requires fastq2")
        if row["layout"] == "SE" and row.get("fastq2", ""):
            raise ValidationError(f"{path}:{row['_line']}: SE run must not provide fastq2")
        if check_fastq:
            for field in ("fastq1", "fastq2"):
                value = row.get(field, "")
                if value and not Path(value).is_file():
                    errors.append(f"{path}:{row['_line']}: {field} does not exist: {value}")
        by_library[row["library_id"]].append(row)

    libraries: list[dict[str, str]] = []
    for library_id, library_rows in sorted(by_library.items()):
        core = {}
        for field in ("study_id", "bio_sample_id", "assay", "condition", "layout"):
            core[field] = consistent(library_rows, field, f"library {library_id}")
        for field in ("batch", "strandedness", "umi_pattern", "ribo_protocol_group", "adapter_3p"):
            values = {row.get(field, "") for row in library_rows if row.get(field, "")}
            if len(values) > 1:
                raise ValidationError(f"library {library_id} has inconsistent {field}: {sorted(values)}")
            core[field] = next(iter(values), "")
        libraries.append(
            {
                "library_id": library_id,
                **core,
                "run_ids": ",".join(row["run_id"] for row in library_rows),
                "n_runs": str(len(library_rows)),
            }
        )

    pair_seen: dict[tuple[str, str, str], str] = {}
    for library in libraries:
        key = (library["study_id"], library["bio_sample_id"], library["assay"])
        if key in pair_seen:
            raise ValidationError(
                f"Biological sample {key[1]} in {key[0]} has multiple logical {key[2]} libraries: "
                f"{pair_seen[key]}, {library['library_id']}"
            )
        pair_seen[key] = library["library_id"]
    return rows, libraries, errors


def validate_contrasts(path: Path, libraries: list[dict[str, str]]) -> tuple[list[dict[str, str]], list[dict[str, object]]]:
    _, rows = read_tsv(path, CONTRAST_REQUIRED)
    conditions: dict[str, set[str]] = defaultdict(set)
    samples: dict[tuple[str, str, str], set[str]] = defaultdict(set)
    pairs: dict[tuple[str, str, str], set[str]] = defaultdict(set)
    for library in libraries:
        study = library["study_id"]
        condition = library["condition"]
        assay = library["assay"]
        bio = library["bio_sample_id"]
        conditions[study].add(condition)
        samples[(study, condition, assay)].add(bio)
        pairs[(study, condition, bio)].add(assay)

    seen: set[tuple[str, str]] = set()
    summaries: list[dict[str, object]] = []
    for row in rows:
        for field in CONTRAST_REQUIRED:
            require_value(row, field, path)
        key = (row["study_id"], row["contrast_id"])
        if not SAFE_ID.fullmatch(row["contrast_id"]):
            raise ValidationError(f"{path}:{row['_line']}: unsafe contrast_id")
        if key in seen:
            raise ValidationError(f"{path}:{row['_line']}: duplicate contrast {key}")
        seen.add(key)
        for flag in ("run_ribo_de", "run_rna_de", "run_dte"):
            if row[flag] not in YES_NO:
                raise ValidationError(f"{path}:{row['_line']}: {flag} must be yes or no")
        if row["test_condition"] == row["reference_condition"]:
            raise ValidationError(f"{path}:{row['_line']}: contrast conditions must differ")
        study_conditions = conditions.get(row["study_id"], set())
        for field in ("test_condition", "reference_condition"):
            if row[field] not in study_conditions:
                raise ValidationError(
                    f"{path}:{row['_line']}: unknown {field} {row[field]!r} for study {row['study_id']}"
                )
        counts: dict[str, int] = {}
        for condition in (row["reference_condition"], row["test_condition"]):
            for assay in ("ribo", "rna"):
                counts[f"{condition}_{assay}"] = len(samples[(row["study_id"], condition, assay)])
            counts[f"{condition}_paired"] = sum(
                pairs[(row["study_id"], condition, bio)] == {"ribo", "rna"}
                for bio in {library["bio_sample_id"] for library in libraries if library["study_id"] == row["study_id"] and library["condition"] == condition}
            )
        if row["run_dte"] == "yes" and min(
            counts[f"{row['reference_condition']}_paired"], counts[f"{row['test_condition']}_paired"]
        ) < 1:
            raise ValidationError(f"{path}:{row['_line']}: DTE contrast has an unpaired condition")
        def replication_level(minimum: int) -> str:
            return "formal" if minimum >= 3 else "exploratory_low_replication" if minimum == 2 else "descriptive"

        summaries.append(
            {
                "study_id": row["study_id"],
                "contrast_id": row["contrast_id"],
                **counts,
                "ribo_analysis_level": replication_level(
                    min(counts[f"{row['reference_condition']}_ribo"], counts[f"{row['test_condition']}_ribo"])
                ),
                "rna_analysis_level": replication_level(
                    min(counts[f"{row['reference_condition']}_rna"], counts[f"{row['test_condition']}_rna"])
                ),
                "dte_analysis_level": replication_level(
                    min(counts[f"{row['reference_condition']}_paired"], counts[f"{row['test_condition']}_paired"])
                ),
            }
        )
    return rows, summaries


def validate_decisions(
    path: Path, libraries: list[dict[str, str]], expected_qc_fingerprint: str
) -> list[dict[str, str]]:
    _, rows = read_tsv(path, DECISION_REQUIRED)
    library_map = {library["library_id"]: library for library in libraries}
    seen: set[str] = set()
    for row in rows:
        library_id = require_value(row, "library_id", path)
        if library_id in seen:
            raise ValidationError(f"{path}:{row['_line']}: duplicate library_id {library_id}")
        seen.add(library_id)
        if library_id not in library_map:
            raise ValidationError(f"{path}:{row['_line']}: unknown library_id {library_id}")
        for flag in ("include_assay_analysis", "include_te"):
            if row[flag] not in YES_NO:
                raise ValidationError(f"{path}:{row['_line']}: {flag} must be yes or no")
        if row["decision_status"] != "reviewed":
            raise ValidationError(f"{path}:{row['_line']}: decision_status must be reviewed")
        if not row["reviewer"]:
            raise ValidationError(f"{path}:{row['_line']}: reviewer is required")
        try:
            date.fromisoformat(row["review_date"])
        except ValueError as exc:
            raise ValidationError(f"{path}:{row['_line']}: review_date must be YYYY-MM-DD") from exc
        if expected_qc_fingerprint and row["qc_fingerprint"] != expected_qc_fingerprint:
            raise ValidationError(f"{path}:{row['_line']}: stale qc_fingerprint")
        library = library_map[library_id]
        if library["assay"] == "ribo" and row["include_assay_analysis"] == "yes":
            lengths = [item for item in row["selected_lengths"].split(",") if item]
            if not lengths:
                raise ValidationError(f"{path}:{row['_line']}: included Ribo library needs selected_lengths")
            if any(not item.isdigit() or not 20 <= int(item) <= 35 for item in lengths):
                raise ValidationError(f"{path}:{row['_line']}: Ribo lengths must be integers in 20..35")
            if row["offset_source"] not in {"sample", "pooled", "manual"}:
                raise ValidationError(f"{path}:{row['_line']}: invalid offset_source")
            offset_pairs: dict[int, int] = {}
            try:
                for pair in row["offsets"].split(","):
                    length, offset = pair.split(":", 1)
                    offset_pairs[int(length)] = int(offset)
            except ValueError as exc:
                raise ValidationError(
                    f"{path}:{row['_line']}: offsets must look like 28:12,29:12"
                ) from exc
            selected = {int(item) for item in lengths}
            if selected - set(offset_pairs):
                raise ValidationError(f"{path}:{row['_line']}: selected lengths are missing offsets")
            if any(not 0 <= value < length for length, value in offset_pairs.items()):
                raise ValidationError(f"{path}:{row['_line']}: offset must be within its read length")
        requested_adapter = library.get("adapter_3p", "")
        if requested_adapter == "infer":
            if row["adapter_review"] != "accepted_inference":
                raise ValidationError(f"{path}:{row['_line']}: inferred adapter requires explicit review")
            if not row["confirmed_adapter_3p"]:
                raise ValidationError(f"{path}:{row['_line']}: confirmed_adapter_3p is required")
        elif requested_adapter in {"", "none", "NA", "."}:
            if row["adapter_review"] not in {"none", "not_applicable"}:
                raise ValidationError(f"{path}:{row['_line']}: adapter_review must be none")
        else:
            if row["adapter_review"] != "explicit":
                raise ValidationError(f"{path}:{row['_line']}: explicit adapter was not acknowledged")
            if row["confirmed_adapter_3p"] != requested_adapter:
                raise ValidationError(f"{path}:{row['_line']}: confirmed adapter differs from samples.tsv")
        if row["include_te"] == "yes" and row["include_assay_analysis"] != "yes":
            raise ValidationError(f"{path}:{row['_line']}: include_te requires include_assay_analysis=yes")
        if row["qc_status"] not in {"pass", "caution", "fail"}:
            raise ValidationError(f"{path}:{row['_line']}: qc_status must be pass, caution, or fail")
        if (
            row["include_assay_analysis"] == "no"
            or row["include_te"] == "no"
            or row["qc_status"] != "pass"
        ) and not row["reason"]:
            raise ValidationError(f"{path}:{row['_line']}: exclusion/caution/fail requires reason")
    missing = sorted(set(library_map) - seen)
    if missing:
        raise ValidationError(f"{path}: decisions missing libraries: {', '.join(missing)}")
    decisions_by_library = {row["library_id"]: row for row in rows}
    paired = {
        (library["study_id"], library["bio_sample_id"], library["assay"]): library["library_id"]
        for library in libraries
    }
    for library in libraries:
        decision = decisions_by_library[library["library_id"]]
        if decision["include_te"] != "yes":
            continue
        other_assay = "rna" if library["assay"] == "ribo" else "ribo"
        counterpart = paired.get((library["study_id"], library["bio_sample_id"], other_assay))
        if not counterpart or decisions_by_library[counterpart]["include_te"] != "yes":
            raise ValidationError(
                f"{path}: include_te must be yes for both assays of {library['bio_sample_id']}"
            )
    return rows


def stable_hash(paths: Iterable[Path]) -> str:
    digest = hashlib.sha256()
    for path in sorted(paths, key=lambda item: str(item)):
        digest.update(path.name.encode())
        digest.update(path.read_bytes())
    return digest.hexdigest()


def write_tsv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", required=True, type=Path)
    parser.add_argument("--contrasts", required=True, type=Path)
    parser.add_argument("--decisions", type=Path)
    parser.add_argument("--expected-qc-fingerprint", default="")
    parser.add_argument("--expected-qc-fingerprint-file", type=Path)
    parser.add_argument("--check-fastq", action="store_true")
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--libraries", required=True, type=Path)
    args = parser.parse_args()

    try:
        _, libraries, path_errors = validate_samples(args.samples, args.check_fastq)
        _, contrast_summaries = validate_contrasts(args.contrasts, libraries)
        decisions = None
        if args.decisions:
            expected = args.expected_qc_fingerprint
            if args.expected_qc_fingerprint_file:
                expected = args.expected_qc_fingerprint_file.read_text(encoding="utf-8").strip()
            decisions = validate_decisions(args.decisions, libraries, expected)
        if path_errors:
            raise ValidationError("\n".join(path_errors))
    except ValidationError as exc:
        print(f"METADATA_VALIDATION_ERROR: {exc}", file=sys.stderr)
        return 2

    fields = [
        "library_id",
        "study_id",
        "bio_sample_id",
        "assay",
        "condition",
        "layout",
        "batch",
        "strandedness",
        "umi_pattern",
        "ribo_protocol_group",
        "adapter_3p",
        "run_ids",
        "n_runs",
    ]
    write_tsv(args.libraries, libraries, fields)
    report = {
        "status": "valid",
        "schema_version": "1.0",
        "input_fingerprint": stable_hash([args.samples, args.contrasts]),
        "n_libraries": len(libraries),
        "n_studies": len({row["study_id"] for row in libraries}),
        "contrasts": contrast_summaries,
        "qc_decisions_reviewed": decisions is not None,
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
