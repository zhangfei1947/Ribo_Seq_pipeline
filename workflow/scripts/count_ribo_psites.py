#!/usr/bin/env python3
"""Count reviewed, length-specific P-sites over gene-level CDS unions."""

from __future__ import annotations

import argparse
from collections import Counter
import csv
from pathlib import Path
import sys

import pysam

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR.parent))

from lib.gtf import iter_gtf  # noqa: E402
from lib.intervals import FeatureInterval, PointIntervalIndex  # noqa: E402
from ribo_qc import biological_strand, psite_position  # noqa: E402


def parse_offsets(value: str) -> dict[int, int]:
    result: dict[int, int] = {}
    for item in value.split(","):
        if not item:
            continue
        length, offset = item.split(":", 1)
        result[int(length)] = int(offset)
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam", required=True, type=Path)
    parser.add_argument("--gene-cds-gtf", required=True, type=Path)
    parser.add_argument("--decisions", required=True, type=Path)
    parser.add_argument("--library-id", required=True)
    parser.add_argument("--strandedness", required=True, choices=("forward", "reverse"))
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--metrics", required=True, type=Path)
    args = parser.parse_args()
    with args.decisions.open(encoding="utf-8") as handle:
        decision_rows = [row for row in csv.DictReader(handle, delimiter="\t") if row["library_id"] == args.library_id]
    if len(decision_rows) != 1:
        raise SystemExit(f"Expected exactly one QC decision for {args.library_id}")
    decision = decision_rows[0]
    selected = {int(item) for item in decision["selected_lengths"].split(",") if item}
    offsets = parse_offsets(decision["offsets"])
    if selected - set(offsets):
        raise SystemExit(f"Missing offsets for lengths: {sorted(selected - set(offsets))}")

    genes: set[str] = set()
    intervals: list[FeatureInterval] = []
    for record in iter_gtf(args.gene_cds_gtf):
        gene_id = record.attributes.get("gene_id", "")
        if gene_id:
            genes.add(gene_id)
            intervals.append(FeatureInterval(record.seqname, record.start, record.end, record.strand, gene_id))
    index = PointIntervalIndex(intervals)
    counts: Counter[str] = Counter()
    fate: Counter[str] = Counter()
    if decision["include_assay_analysis"] == "yes":
        with pysam.AlignmentFile(args.bam, "rb") as bam:
            for read in bam.fetch(until_eof=True):
                length = read.query_length or 0
                if length not in selected:
                    fate["length_excluded"] += 1
                    continue
                position = psite_position(read, offsets[length])
                if position is None:
                    fate["offset_unmappable"] += 1
                    continue
                strand = biological_strand(read, args.strandedness)
                seqname = bam.get_reference_name(read.reference_id)
                overlaps = {feature.feature_id for feature in index.query(seqname, position, strand)}
                if len(overlaps) == 1:
                    counts[next(iter(overlaps))] += 1
                    fate["assigned"] += 1
                elif overlaps:
                    fate["cross_gene_ambiguous"] += 1
                else:
                    fate["outside_cds"] += 1
    else:
        fate["excluded_by_review"] = 1
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["gene_id", "count"])
        for gene_id in sorted(genes):
            writer.writerow([gene_id, counts[gene_id]])
    args.metrics.parent.mkdir(parents=True, exist_ok=True)
    with args.metrics.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["fate", "count"])
        for name, count in sorted(fate.items()):
            writer.writerow([name, count])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
