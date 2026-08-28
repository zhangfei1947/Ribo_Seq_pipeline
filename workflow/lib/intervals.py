"""Small strand-aware point interval indexes for GTF-sized annotations."""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from typing import Iterable


@dataclass(frozen=True)
class FeatureInterval:
    seqname: str
    start: int
    end: int
    strand: str
    feature_id: str
    payload: object = None


class PointIntervalIndex:
    """Bin 1-based closed intervals for fast point overlap queries."""

    def __init__(self, intervals: Iterable[FeatureInterval], bin_size: int = 16_384):
        self.bin_size = bin_size
        self._bins: dict[tuple[str, str, int], list[FeatureInterval]] = defaultdict(list)
        for interval in intervals:
            if interval.start < 1 or interval.end < interval.start:
                raise ValueError(f"Invalid interval: {interval}")
            first = interval.start // bin_size
            last = interval.end // bin_size
            for bin_number in range(first, last + 1):
                self._bins[(interval.seqname, interval.strand, bin_number)].append(interval)

    def query(self, seqname: str, position: int, strand: str) -> list[FeatureInterval]:
        return [
            interval
            for interval in self._bins.get((seqname, strand, position // self.bin_size), ())
            if interval.start <= position <= interval.end
        ]

