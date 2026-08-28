from __future__ import annotations

import csv
from pathlib import Path
import tempfile
import unittest

from workflow.scripts.validate_metadata import ValidationError, validate_contrasts, validate_samples


SAMPLE_FIELDS = [
    "study_id", "run_id", "library_id", "bio_sample_id", "assay", "condition",
    "fastq1", "fastq2", "layout", "batch", "strandedness", "umi_pattern",
    "ribo_protocol_group"
]


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


class MetadataTests(unittest.TestCase):
    def sample_rows(self):
        rows = []
        for condition in ("control", "treated"):
            for index in (1, 2):
                bio = f"{condition}_{index}"
                rows.extend([
                    {"study_id": "s1", "run_id": f"ribo_{bio}", "library_id": f"ribo_{bio}", "bio_sample_id": bio, "assay": "ribo", "condition": condition, "fastq1": "/x", "fastq2": "", "layout": "SE", "batch": "b1", "strandedness": "forward", "umi_pattern": "none", "ribo_protocol_group": "p1"},
                    {"study_id": "s1", "run_id": f"rna_{bio}", "library_id": f"rna_{bio}", "bio_sample_id": bio, "assay": "rna", "condition": condition, "fastq1": "/x", "fastq2": "/y", "layout": "PE", "batch": "b1", "strandedness": "reverse", "umi_pattern": "none", "ribo_protocol_group": ""},
                ])
        return rows

    def test_valid_paired_design_is_exploratory_with_two_replicates(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            samples = tmp / "samples.tsv"
            contrasts = tmp / "contrasts.tsv"
            write_tsv(samples, SAMPLE_FIELDS, self.sample_rows())
            write_tsv(contrasts, ["study_id", "contrast_id", "test_condition", "reference_condition", "run_ribo_de", "run_rna_de", "run_dte"], [{"study_id": "s1", "contrast_id": "t_vs_c", "test_condition": "treated", "reference_condition": "control", "run_ribo_de": "yes", "run_rna_de": "yes", "run_dte": "yes"}])
            _, libraries, errors = validate_samples(samples, False)
            self.assertEqual(errors, [])
            _, summaries = validate_contrasts(contrasts, libraries)
            self.assertEqual(summaries[0]["dte_analysis_level"], "exploratory_low_replication")
            self.assertEqual(summaries[0]["ribo_analysis_level"], "exploratory_low_replication")

    def test_ribo_pe_is_rejected(self):
        rows = self.sample_rows()
        rows[0]["layout"] = "PE"
        rows[0]["fastq2"] = "/y"
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "samples.tsv"
            write_tsv(path, SAMPLE_FIELDS, rows)
            with self.assertRaises(ValidationError):
                validate_samples(path, False)


if __name__ == "__main__":
    unittest.main()
