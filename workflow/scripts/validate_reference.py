#!/usr/bin/env python3
"""Fail fast when the FlyBase GTF and prebuilt alignment indexes do not match policy."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR.parent))

from lib.gtf import iter_gtf  # noqa: E402


def star_parameters(directory: Path) -> dict[str, str]:
    required = ("Genome", "SA", "SAindex", "genomeParameters.txt", "chrName.txt")
    missing = [name for name in required if not (directory / name).is_file()]
    if missing:
        raise ValueError(f"{directory}: missing STAR files: {', '.join(missing)}")
    parameters: dict[str, str] = {}
    for line in (directory / "genomeParameters.txt").read_text(encoding="utf-8").splitlines():
        fields = line.split()
        if len(fields) >= 2:
            parameters[fields[0]] = " ".join(fields[1:])
    return parameters


def chromosomes(directory: Path) -> set[str]:
    return {
        line.strip()
        for line in (directory / "chrName.txt").read_text(encoding="utf-8").splitlines()
        if line.strip()
    }


def bowtie_files(prefix: Path) -> list[Path]:
    suffixes = ("1", "2", "3", "4", "rev.1", "rev.2")
    small = [Path(f"{prefix}.{suffix}.bt2") for suffix in suffixes]
    large = [Path(f"{prefix}.{suffix}.bt2l") for suffix in suffixes]
    if all(path.is_file() for path in small):
        return small
    if all(path.is_file() for path in large):
        return large
    raise ValueError(f"Incomplete Bowtie2 index for prefix {prefix}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gtf", required=True, type=Path)
    parser.add_argument("--ribo-star", required=True, type=Path)
    parser.add_argument("--rna-star", required=True, type=Path)
    parser.add_argument("--contaminant-prefix", required=True, type=Path)
    parser.add_argument("--ribo-overhang", required=True, type=int)
    parser.add_argument("--rna-overhang", required=True, type=int)
    parser.add_argument("--sa-index-bases", default=12, type=int)
    parser.add_argument("--release", required=True)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    ribo_params = star_parameters(args.ribo_star)
    rna_params = star_parameters(args.rna_star)
    checks = (
        ("Ribo STAR sjdbOverhang", ribo_params.get("sjdbOverhang"), str(args.ribo_overhang)),
        ("RNA STAR sjdbOverhang", rna_params.get("sjdbOverhang"), str(args.rna_overhang)),
        ("Ribo STAR genomeSAindexNbases", ribo_params.get("genomeSAindexNbases"), str(args.sa_index_bases)),
        ("RNA STAR genomeSAindexNbases", rna_params.get("genomeSAindexNbases"), str(args.sa_index_bases)),
    )
    for label, observed, expected in checks:
        if observed != expected:
            raise ValueError(f"{label}: expected {expected}, observed {observed!r}")
    contaminant_files = bowtie_files(args.contaminant_prefix)

    gtf_chromosomes: set[str] = set()
    features: set[str] = set()
    gene_symbols = 0
    records = 0
    for record in iter_gtf(args.gtf):
        records += 1
        gtf_chromosomes.add(record.seqname)
        features.add(record.feature)
        gene_symbols += bool(record.attributes.get("gene_symbol"))
    if not records or not {"gene", "exon", "CDS"}.issubset(features):
        raise ValueError("GTF must contain gene, exon, and CDS records")
    if not gene_symbols:
        raise ValueError("GTF contains no gene_symbol attributes")
    for label, directory in (("Ribo STAR", args.ribo_star), ("RNA STAR", args.rna_star)):
        missing = sorted(gtf_chromosomes - chromosomes(directory))
        if missing:
            raise ValueError(f"{label} is missing GTF chromosomes: {', '.join(missing[:10])}")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(
            {
                "status": "valid",
                "release": args.release,
                "gtf": str(args.gtf),
                "gtf_sha256": sha256(args.gtf),
                "gtf_records": records,
                "gtf_chromosomes": sorted(gtf_chromosomes),
                "ribo_star_parameters": ribo_params,
                "rna_star_parameters": rna_params,
                "contaminant_index_files": [str(path) for path in contaminant_files],
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

