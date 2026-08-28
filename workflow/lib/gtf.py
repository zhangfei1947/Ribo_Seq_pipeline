"""Small, dependency-free GTF parsing and interval utilities."""

from __future__ import annotations

from dataclasses import dataclass
import gzip
from pathlib import Path
import re
from typing import Iterable, Iterator, Mapping, TextIO


ATTRIBUTE_RE = re.compile(r"\s*([^\s;]+)\s+(?:\"([^\"]*)\"|([^;\s]+))\s*;?")


class GTFError(ValueError):
    """Raised for malformed or internally inconsistent GTF input."""


@dataclass(frozen=True)
class GTFRecord:
    seqname: str
    source: str
    feature: str
    start: int
    end: int
    score: str
    strand: str
    frame: str
    attributes: Mapping[str, str]


def open_text(path: str | Path) -> TextIO:
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("rt", encoding="utf-8")


def parse_attributes(raw: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    pos = 0
    while pos < len(raw):
        if raw[pos:].strip() == "":
            break
        match = ATTRIBUTE_RE.match(raw, pos)
        if not match:
            fragment = raw[pos : pos + 80]
            raise GTFError(f"Cannot parse GTF attribute fragment: {fragment!r}")
        key = match.group(1)
        value = match.group(2) if match.group(2) is not None else match.group(3)
        if key in attrs and attrs[key] != value:
            # Repeated tags are legal. Preserve them without losing information.
            attrs[key] = f"{attrs[key]},{value}"
        else:
            attrs[key] = value or ""
        pos = match.end()
    return attrs


def iter_gtf(path: str | Path) -> Iterator[GTFRecord]:
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise GTFError(f"{path}:{line_number}: expected 9 columns, found {len(fields)}")
            try:
                start = int(fields[3])
                end = int(fields[4])
            except ValueError as exc:
                raise GTFError(f"{path}:{line_number}: invalid coordinates") from exc
            if start < 1 or end < start:
                raise GTFError(f"{path}:{line_number}: invalid 1-based closed interval {start}-{end}")
            yield GTFRecord(
                seqname=fields[0],
                source=fields[1],
                feature=fields[2],
                start=start,
                end=end,
                score=fields[5],
                strand=fields[6],
                frame=fields[7],
                attributes=parse_attributes(fields[8]),
            )


def merge_intervals(intervals: Iterable[tuple[int, int]]) -> list[tuple[int, int]]:
    """Merge 1-based closed intervals, including directly adjacent intervals."""
    ordered = sorted(intervals)
    if not ordered:
        return []
    merged: list[list[int]] = [[ordered[0][0], ordered[0][1]]]
    for start, end in ordered[1:]:
        current = merged[-1]
        if start <= current[1] + 1:
            current[1] = max(current[1], end)
        else:
            merged.append([start, end])
    return [(start, end) for start, end in merged]


def union_length(intervals: Iterable[tuple[int, int]]) -> int:
    return sum(end - start + 1 for start, end in merge_intervals(intervals))

