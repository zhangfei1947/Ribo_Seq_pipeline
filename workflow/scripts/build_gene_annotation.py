#!/usr/bin/env python3
"""Build deterministic FlyBase gene and representative-transcript annotations."""

from __future__ import annotations

import argparse
from collections import defaultdict
import csv
from dataclasses import dataclass, field
from pathlib import Path
import sys

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR.parent))

from lib.gtf import GTFError, iter_gtf, merge_intervals, union_length  # noqa: E402


GENE_FIELDS = [
    "gene_id",
    "gene_symbol",
    "seqname",
    "start",
    "end",
    "strand",
    "gene_span_length",
    "exon_union_length",
    "cds_union_length",
    "transcript_count",
    "has_cds",
    "gene_biotype",
    "annotation_release",
]

REP_FIELDS = [
    "gene_id",
    "gene_symbol",
    "transcript_id",
    "transcript_symbol",
    "seqname",
    "start",
    "end",
    "strand",
    "transcript_span_length",
    "cds_union_length",
    "selection_reason",
    "annotation_release",
]


@dataclass
class GeneInfo:
    gene_id: str
    symbol: str = ""
    seqname: str = ""
    start: int | None = None
    end: int | None = None
    strand: str = ""
    strands: set[str] = field(default_factory=set)
    biotype: str = ""
    exons: list[tuple[int, int, str]] = field(default_factory=list)
    cds: list[tuple[int, int, str]] = field(default_factory=list)
    transcripts: set[str] = field(default_factory=set)


@dataclass
class TranscriptInfo:
    transcript_id: str
    gene_id: str
    symbol: str = ""
    transcript_symbol: str = ""
    seqname: str = ""
    start: int | None = None
    end: int | None = None
    strand: str = ""
    strands: set[str] = field(default_factory=set)
    cds: list[tuple[int, int, str]] = field(default_factory=list)
    canonical: bool = False


def first_attr(attrs: dict[str, str] | object, *keys: str) -> str:
    for key in keys:
        value = attrs.get(key, "")  # type: ignore[union-attr]
        if value:
            return value
    return ""


def update_location(obj: GeneInfo | TranscriptInfo, seqname: str, start: int, end: int, strand: str) -> None:
    if obj.seqname and obj.seqname != seqname:
        raise GTFError(f"{getattr(obj, 'gene_id', getattr(obj, 'transcript_id', '?'))} spans chromosomes")
    obj.seqname = seqname
    if strand in {"+", "-"}:
        obj.strands.add(strand)
    if len(obj.strands) == 1:
        obj.strand = next(iter(obj.strands))
    elif len(obj.strands) > 1:
        # FlyBase contains genuine trans-spliced loci (for example mod(mdg4))
        # whose records span both strands. Preserve that fact at gene/transcript
        # level while retaining the strand of each exon/CDS interval below.
        obj.strand = "."
    elif not obj.strand:
        obj.strand = strand
    obj.start = start if obj.start is None else min(obj.start, start)
    obj.end = end if obj.end is None else max(obj.end, end)


def set_consistent(current: str, new: str, label: str) -> str:
    if current and new and current != new:
        raise GTFError(f"Conflicting {label}: {current!r} vs {new!r}")
    return current or new


def build(gtf: Path) -> tuple[dict[str, GeneInfo], dict[str, TranscriptInfo]]:
    genes: dict[str, GeneInfo] = {}
    transcripts: dict[str, TranscriptInfo] = {}

    for record in iter_gtf(gtf):
        attrs = dict(record.attributes)
        gene_id = attrs.get("gene_id", "")
        if not gene_id:
            continue
        gene = genes.setdefault(gene_id, GeneInfo(gene_id=gene_id))
        symbol = first_attr(attrs, "gene_symbol", "gene_name")
        gene.symbol = set_consistent(gene.symbol, symbol, f"gene_symbol for {gene_id}")
        gene.biotype = set_consistent(
            gene.biotype,
            first_attr(attrs, "gene_biotype", "gene_type", "biotype"),
            f"gene biotype for {gene_id}",
        )
        update_location(gene, record.seqname, record.start, record.end, record.strand)

        transcript_id = attrs.get("transcript_id", "")
        if transcript_id:
            gene.transcripts.add(transcript_id)
            transcript = transcripts.setdefault(
                transcript_id, TranscriptInfo(transcript_id=transcript_id, gene_id=gene_id)
            )
            if transcript.gene_id != gene_id:
                raise GTFError(f"Transcript {transcript_id} maps to multiple genes")
            transcript.symbol = set_consistent(
                transcript.symbol, symbol, f"gene_symbol for transcript {transcript_id}"
            )
            transcript.transcript_symbol = set_consistent(
                transcript.transcript_symbol,
                first_attr(attrs, "transcript_symbol", "transcript_name"),
                f"transcript_symbol for {transcript_id}",
            )
            update_location(transcript, record.seqname, record.start, record.end, record.strand)
            tags = ",".join(
                first_attr(attrs, "tag", "transcript_status", "canonical").lower().split(",")
            )
            if "canonical" in tags or attrs.get("canonical", "").lower() in {"1", "true", "yes"}:
                transcript.canonical = True

        if record.feature == "exon":
            gene.exons.append((record.start, record.end, record.strand))
        elif record.feature == "CDS":
            gene.cds.append((record.start, record.end, record.strand))
            if transcript_id:
                transcripts[transcript_id].cds.append((record.start, record.end, record.frame))

    if not genes:
        raise GTFError(f"No gene_id records found in {gtf}")
    return genes, transcripts


def choose_representatives(
    genes: dict[str, GeneInfo], transcripts: dict[str, TranscriptInfo]
) -> dict[str, tuple[TranscriptInfo, str]]:
    selected: dict[str, tuple[TranscriptInfo, str]] = {}
    for gene_id, gene in genes.items():
        candidates = [
            transcripts[tid]
            for tid in gene.transcripts
            if tid in transcripts and transcripts[tid].cds and transcripts[tid].strand in {"+", "-"}
        ]
        if not candidates:
            continue
        canonical = [tx for tx in candidates if tx.canonical]
        pool = canonical or candidates
        pool.sort(
            key=lambda tx: (
                -union_length((start, end) for start, end, _ in tx.cds),
                -((tx.end or 0) - (tx.start or 1) + 1),
                tx.transcript_id,
            )
        )
        selected[gene_id] = (pool[0], "canonical" if canonical else "longest_cds")
    return selected


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def gtf_attributes(gene_id: str, gene_symbol: str) -> str:
    escaped_id = gene_id.replace('"', '')
    escaped_symbol = gene_symbol.replace('"', '')
    return f'gene_id "{escaped_id}"; gene_symbol "{escaped_symbol}";'


def write_union_gtf(
    path: Path, genes: dict[str, GeneInfo], feature: str, interval_attribute: str
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for gene_id in sorted(genes):
            gene = genes[gene_id]
            intervals_by_strand: dict[str, list[tuple[int, int]]] = defaultdict(list)
            for start, end, strand in getattr(gene, interval_attribute):
                intervals_by_strand[strand].append((start, end))
            for strand in sorted(intervals_by_strand):
                for start, end in merge_intervals(intervals_by_strand[strand]):
                    fields = [
                        gene.seqname,
                        "RiboSeqPipeline",
                        feature,
                        str(start),
                        str(end),
                        ".",
                        strand,
                        ".",
                        gtf_attributes(gene_id, gene.symbol),
                    ]
                    handle.write("\t".join(fields) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gtf", required=True, type=Path)
    parser.add_argument("--annotation-release", required=True)
    parser.add_argument("--all-genes", required=True, type=Path)
    parser.add_argument("--cds-genes", required=True, type=Path)
    parser.add_argument("--representative-transcripts", required=True, type=Path)
    parser.add_argument("--gene-cds-gtf", required=True, type=Path)
    parser.add_argument("--gene-exon-gtf", required=True, type=Path)
    parser.add_argument("--representative-cds", required=True, type=Path)
    args = parser.parse_args()

    genes, transcripts = build(args.gtf)
    representatives = choose_representatives(genes, transcripts)
    rows: list[dict[str, object]] = []
    for gene_id in sorted(genes):
        gene = genes[gene_id]
        if gene.start is None or gene.end is None:
            raise GTFError(f"Gene {gene_id} has no coordinates")
        rows.append(
            {
                "gene_id": gene_id,
                "gene_symbol": gene.symbol,
                "seqname": gene.seqname,
                "start": gene.start,
                "end": gene.end,
                "strand": gene.strand,
                "gene_span_length": gene.end - gene.start + 1,
                "exon_union_length": union_length((start, end) for start, end, _ in gene.exons),
                "cds_union_length": union_length((start, end) for start, end, _ in gene.cds),
                "transcript_count": len(gene.transcripts),
                "has_cds": "yes" if gene.cds else "no",
                "gene_biotype": gene.biotype,
                "annotation_release": args.annotation_release,
            }
        )
    write_tsv(args.all_genes, GENE_FIELDS, rows)
    write_tsv(args.cds_genes, GENE_FIELDS, [row for row in rows if row["has_cds"] == "yes"])

    rep_rows: list[dict[str, object]] = []
    for gene_id in sorted(representatives):
        tx, reason = representatives[gene_id]
        gene = genes[gene_id]
        if tx.start is None or tx.end is None:
            continue
        rep_rows.append(
            {
                "gene_id": gene_id,
                "gene_symbol": gene.symbol,
                "transcript_id": tx.transcript_id,
                "transcript_symbol": tx.transcript_symbol,
                "seqname": tx.seqname,
                "start": tx.start,
                "end": tx.end,
                "strand": tx.strand,
                "transcript_span_length": tx.end - tx.start + 1,
                "cds_union_length": union_length((start, end) for start, end, _ in tx.cds),
                "selection_reason": reason,
                "annotation_release": args.annotation_release,
            }
        )
    write_tsv(args.representative_transcripts, REP_FIELDS, rep_rows)
    write_union_gtf(args.gene_cds_gtf, genes, "CDS", "cds")
    write_union_gtf(args.gene_exon_gtf, genes, "exon", "exons")

    block_rows: list[dict[str, object]] = []
    for gene_id in sorted(representatives):
        tx, _ = representatives[gene_id]
        ordered = sorted(tx.cds, key=lambda item: item[0], reverse=tx.strand == "-")
        for block_index, (start, end, phase) in enumerate(ordered, start=1):
            block_rows.append(
                {
                    "gene_id": gene_id,
                    "gene_symbol": genes[gene_id].symbol,
                    "transcript_id": tx.transcript_id,
                    "block_index": block_index,
                    "seqname": tx.seqname,
                    "start": start,
                    "end": end,
                    "strand": tx.strand,
                    "phase": phase,
                }
            )
    write_tsv(
        args.representative_cds,
        [
            "gene_id",
            "gene_symbol",
            "transcript_id",
            "block_index",
            "seqname",
            "start",
            "end",
            "strand",
            "phase",
        ],
        block_rows,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
