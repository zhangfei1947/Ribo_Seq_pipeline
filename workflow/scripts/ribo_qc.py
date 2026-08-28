#!/usr/bin/env python3
"""Infer length-specific P-site offsets and produce Ribo-seq QC tables."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import json
from pathlib import Path
import sys

import pysam

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR.parent))

from lib.intervals import FeatureInterval, PointIntervalIndex  # noqa: E402
from lib.gtf import iter_gtf  # noqa: E402


def parse_range(value: str) -> list[int]:
    if "-" in value:
        left, right = value.split("-", 1)
        return list(range(int(left), int(right) + 1))
    return sorted({int(item) for item in value.split(",") if item})


def load_representative_cds(path: Path) -> tuple[PointIntervalIndex, dict[str, int]]:
    by_tx: dict[str, list[dict[str, str]]] = defaultdict(list)
    with path.open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            by_tx[row["transcript_id"]].append(row)
    intervals: list[FeatureInterval] = []
    transcript_lengths: dict[str, int] = {}
    for transcript_id, rows in by_tx.items():
        rows.sort(key=lambda row: int(row["block_index"]))
        first_phase = int(rows[0]["phase"]) if rows[0]["phase"] in {"0", "1", "2"} else 0
        cumulative = 0
        for row in rows:
            start, end = int(row["start"]), int(row["end"])
            intervals.append(
                FeatureInterval(
                    row["seqname"],
                    start,
                    end,
                    row["strand"],
                    row["gene_id"],
                    {
                        "transcript_id": transcript_id,
                        "cumulative": cumulative,
                        "strand": row["strand"],
                        "start": start,
                        "end": end,
                        "first_phase": first_phase,
                    },
                )
            )
            cumulative += end - start + 1
        transcript_lengths[transcript_id] = cumulative
    return PointIntervalIndex(intervals), transcript_lengths


def load_gene_features(path: Path) -> tuple[PointIntervalIndex, set[str]]:
    intervals: list[FeatureInterval] = []
    genes: set[str] = set()
    for record in iter_gtf(path):
        gene_id = record.attributes.get("gene_id", "")
        if gene_id:
            genes.add(gene_id)
            intervals.append(FeatureInterval(record.seqname, record.start, record.end, record.strand, gene_id))
    return PointIntervalIndex(intervals), genes


def psite_position(read: pysam.AlignedSegment, offset: int) -> int | None:
    positions = read.get_reference_positions(full_length=True)
    index = len(positions) - offset - 1 if read.is_reverse else offset
    if index < 0 or index >= len(positions):
        return None
    position = positions[index]
    return None if position is None else position + 1


def biological_strand(read: pysam.AlignedSegment, strandedness: str) -> str:
    alignment_strand = "-" if read.is_reverse else "+"
    if strandedness == "forward":
        return alignment_strand
    return "+" if alignment_strand == "-" else "-"


def transcript_coordinate(feature: FeatureInterval, position: int) -> tuple[str, int, int]:
    payload = feature.payload
    if not isinstance(payload, dict):
        raise TypeError("Representative CDS interval is missing payload")
    within = position - feature.start if payload["strand"] == "+" else feature.end - position
    return str(payload["transcript_id"]), int(payload["cumulative"]) + within, int(payload["first_phase"])


def unique_cds_coordinate(
    index: PointIntervalIndex, lengths: dict[str, int], seqname: str, position: int, strand: str
) -> tuple[int, int, int] | None:
    candidates = index.query(seqname, position, strand)
    mapped = {(feature.feature_id, *transcript_coordinate(feature, position)) for feature in candidates}
    if len(mapped) != 1:
        return None
    _, transcript_id, coordinate, first_phase = next(iter(mapped))
    frame = (int(coordinate) - int(first_phase)) % 3
    return int(coordinate), lengths[str(transcript_id)], frame


def write_rows(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam", required=True, type=Path)
    parser.add_argument("--representative-cds", required=True, type=Path)
    parser.add_argument("--gene-cds-gtf", required=True, type=Path)
    parser.add_argument("--gene-exon-gtf", required=True, type=Path)
    parser.add_argument("--candidate-lengths", default="20-35")
    parser.add_argument("--offset-min", default=6, type=int)
    parser.add_argument("--offset-max", default=20, type=int)
    parser.add_argument("--min-reads-per-length", default=100, type=int)
    parser.add_argument("--strandedness", required=True, choices=("forward", "reverse"))
    parser.add_argument("--length-metrics", required=True, type=Path)
    parser.add_argument("--offsets", required=True, type=Path)
    parser.add_argument("--offset-candidates", required=True, type=Path)
    parser.add_argument("--metagene", required=True, type=Path)
    parser.add_argument("--provisional-counts", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    args = parser.parse_args()

    allowed_lengths = set(parse_range(args.candidate_lengths))
    rep_index, transcript_lengths = load_representative_cds(args.representative_cds)
    gene_cds_index, coding_genes = load_gene_features(args.gene_cds_gtf)
    gene_exon_index, _ = load_gene_features(args.gene_exon_gtf)
    length_totals: Counter[int] = Counter()
    candidate_frames: dict[tuple[int, int], Counter[int]] = defaultdict(Counter)

    with pysam.AlignmentFile(args.bam, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            length = read.query_length or 0
            length_totals[length] += 1
            if length not in allowed_lengths:
                continue
            strand = biological_strand(read, args.strandedness)
            seqname = bam.get_reference_name(read.reference_id)
            for offset in range(args.offset_min, args.offset_max + 1):
                position = psite_position(read, offset)
                if position is None:
                    continue
                coordinate = unique_cds_coordinate(rep_index, transcript_lengths, seqname, position, strand)
                if coordinate is not None:
                    candidate_frames[(length, offset)][coordinate[2]] += 1

    offset_rows: list[dict[str, object]] = []
    candidate_rows: list[dict[str, object]] = []
    chosen: dict[int, int] = {}
    for length in sorted(allowed_lengths):
        scores: list[tuple[float, int, int, Counter[int]]] = []
        for offset in range(args.offset_min, args.offset_max + 1):
            frames = candidate_frames[(length, offset)]
            total = sum(frames.values())
            score = frames[0] / total if total else 0.0
            scores.append((score, total, -offset, frames))
            candidate_rows.append(
                {
                    "read_length": length,
                    "offset": offset,
                    "n_representative_cds": total,
                    "frame0": frames[0],
                    "frame1": frames[1],
                    "frame2": frames[2],
                    "frame0_fraction": f"{score:.6f}",
                }
            )
        scores.sort(reverse=True)
        best_score, best_n, negative_offset, best_frames = scores[0]
        second_score = scores[1][0]
        status = "recommended" if best_n >= args.min_reads_per_length else "insufficient_reads"
        offset = -negative_offset
        if status == "recommended":
            chosen[length] = offset
        offset_rows.append(
            {
                "read_length": length,
                "recommended_offset": offset,
                "n_representative_cds": best_n,
                "frame0": best_frames[0],
                "frame1": best_frames[1],
                "frame2": best_frames[2],
                "frame0_fraction": f"{best_score:.6f}",
                "second_best_margin": f"{best_score - second_score:.6f}",
                "status": status,
            }
        )

    metagene_counts: Counter[tuple[int, str, int]] = Counter()
    selected_frames: dict[int, Counter[int]] = defaultdict(Counter)
    region_counts: Counter[str] = Counter()
    provisional_gene_counts: Counter[str] = Counter()
    with pysam.AlignmentFile(args.bam, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            length = read.query_length or 0
            if length not in chosen:
                continue
            position = psite_position(read, chosen[length])
            if position is None:
                continue
            strand = biological_strand(read, args.strandedness)
            seqname = bam.get_reference_name(read.reference_id)
            cds_genes = {feature.feature_id for feature in gene_cds_index.query(seqname, position, strand)}
            exon_genes = {feature.feature_id for feature in gene_exon_index.query(seqname, position, strand)}
            if cds_genes:
                region_counts["cds"] += 1
                if len(cds_genes) == 1:
                    provisional_gene_counts[next(iter(cds_genes))] += 1
                else:
                    region_counts["cross_gene_ambiguous"] += 1
            elif exon_genes:
                region_counts["exonic_non_cds"] += 1
            else:
                region_counts["non_exonic"] += 1
            coordinate = unique_cds_coordinate(rep_index, transcript_lengths, seqname, position, strand)
            if coordinate is None:
                continue
            cds_coordinate, cds_length, frame = coordinate
            selected_frames[length][frame] += 1
            relative_start = cds_coordinate
            relative_stop = cds_coordinate - cds_length
            if -60 <= relative_start <= 90:
                metagene_counts[(length, "start", relative_start)] += 1
            if -90 <= relative_stop <= 60:
                metagene_counts[(length, "stop", relative_stop)] += 1

    length_rows = []
    all_reads = sum(length_totals.values())
    for length, count in sorted(length_totals.items()):
        frames = selected_frames[length]
        frame_total = sum(frames.values())
        length_rows.append(
            {
                "read_length": length,
                "unique_alignments": count,
                "fraction_unique": f"{count / all_reads:.6f}" if all_reads else "0",
                "selected_offset": chosen.get(length, ""),
                "representative_cds_psites": frame_total,
                "frame0_fraction": f"{frames[0] / frame_total:.6f}" if frame_total else "",
            }
        )
    metagene_rows = [
        {"read_length": length, "landmark": landmark, "position": position, "count": count}
        for (length, landmark, position), count in sorted(metagene_counts.items())
    ]
    write_rows(
        args.length_metrics,
        ["read_length", "unique_alignments", "fraction_unique", "selected_offset", "representative_cds_psites", "frame0_fraction"],
        length_rows,
    )
    write_rows(
        args.offsets,
        ["read_length", "recommended_offset", "n_representative_cds", "frame0", "frame1", "frame2", "frame0_fraction", "second_best_margin", "status"],
        offset_rows,
    )
    write_rows(
        args.offset_candidates,
        ["read_length", "offset", "n_representative_cds", "frame0", "frame1", "frame2", "frame0_fraction"],
        candidate_rows,
    )
    write_rows(args.metagene, ["read_length", "landmark", "position", "count"], metagene_rows)
    write_rows(
        args.provisional_counts,
        ["gene_id", "count"],
        [{"gene_id": gene_id, "count": provisional_gene_counts[gene_id]} for gene_id in sorted(coding_genes)],
    )
    region_total = region_counts["cds"] + region_counts["exonic_non_cds"] + region_counts["non_exonic"]
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    args.summary.write_text(
        json.dumps(
            {
                "status": "complete",
                "unique_alignments": all_reads,
                "candidate_lengths": sorted(allowed_lengths),
                "recommended_offsets": {str(length): offset for length, offset in chosen.items()},
                "psite_regions": dict(region_counts),
                "cds_fraction": region_counts["cds"] / region_total if region_total else None,
                "warning": "Offsets are recommendations and must be frozen in reviewed qc_decisions.tsv",
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
