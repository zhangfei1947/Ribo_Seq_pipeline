from __future__ import annotations

import csv
from pathlib import Path
import subprocess
import tempfile
import unittest

from workflow.lib.gtf import merge_intervals, parse_attributes, union_length


ROOT = Path(__file__).resolve().parents[1]


class GTFTests(unittest.TestCase):
    def test_attributes(self):
        attrs = parse_attributes('gene_id "FBgn0031081"; gene_symbol "Nep3";')
        self.assertEqual(attrs["gene_id"], "FBgn0031081")
        self.assertEqual(attrs["gene_symbol"], "Nep3")

    def test_closed_interval_union(self):
        self.assertEqual(merge_intervals([(1, 3), (4, 5), (10, 10)]), [(1, 5), (10, 10)])
        self.assertEqual(union_length([(1, 3), (3, 5)]), 5)

    def test_annotation_export(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            all_genes = tmp / "all.tsv"
            cds_genes = tmp / "cds.tsv"
            reps = tmp / "reps.tsv"
            cds_gtf = tmp / "cds.gtf"
            exon_gtf = tmp / "exon.gtf"
            rep_cds = tmp / "rep_cds.tsv"
            subprocess.run(
                [
                    "python3",
                    str(ROOT / "workflow/scripts/build_gene_annotation.py"),
                    "--gtf",
                    str(ROOT / "tests/fixtures/mini.gtf"),
                    "--annotation-release",
                    "test_release",
                    "--all-genes",
                    str(all_genes),
                    "--cds-genes",
                    str(cds_genes),
                    "--representative-transcripts",
                    str(reps),
                    "--gene-cds-gtf",
                    str(cds_gtf),
                    "--gene-exon-gtf",
                    str(exon_gtf),
                    "--representative-cds",
                    str(rep_cds),
                ],
                check=True,
            )
            with all_genes.open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(rows[0]["gene_symbol"], "GeneA")
            self.assertEqual(rows[0]["gene_span_length"], "101")
            self.assertEqual(rows[0]["exon_union_length"], "101")
            # CDS intervals 110-130, 150-180, and 120-170 merge to 110-180.
            self.assertEqual(rows[0]["cds_union_length"], "71")
            self.assertEqual(rows[1]["has_cds"], "no")
            with reps.open() as handle:
                rep_rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(rep_rows[0]["transcript_id"], "FBtr0000001")
            self.assertIn('gene_symbol "GeneA"', cds_gtf.read_text())
            with rep_cds.open() as handle:
                blocks = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual([row["block_index"] for row in blocks], ["1", "2"])

    def test_multistrand_gene_preserves_interval_strands(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            gtf = tmp / "multistrand.gtf"
            gtf.write_text(
                "\n".join(
                    [
                        '3R\tFlyBase\tgene\t100\t400\t.\t-\t.\tgene_id "FBgnMixed"; gene_symbol "Mixed";',
                        '3R\tFlyBase\tmRNA\t100\t200\t.\t-\t.\tgene_id "FBgnMixed"; gene_symbol "Mixed"; transcript_id "FBtrMinus";',
                        '3R\tFlyBase\texon\t100\t200\t.\t-\t.\tgene_id "FBgnMixed"; gene_symbol "Mixed"; transcript_id "FBtrMinus";',
                        '3R\tFlyBase\tCDS\t120\t180\t.\t-\t0\tgene_id "FBgnMixed"; gene_symbol "Mixed"; transcript_id "FBtrMinus";',
                        '3R\tFlyBase\tmRNA\t250\t400\t.\t.\t.\tgene_id "FBgnMixed"; gene_symbol "Mixed"; transcript_id "FBtrMixed";',
                        '3R\tFlyBase\texon\t250\t300\t.\t+\t.\tgene_id "FBgnMixed"; gene_symbol "Mixed"; transcript_id "FBtrMixed";',
                        '3R\tFlyBase\tCDS\t260\t300\t.\t+\t0\tgene_id "FBgnMixed"; gene_symbol "Mixed"; transcript_id "FBtrMixed";',
                        '3R\tFlyBase\texon\t350\t400\t.\t-\t.\tgene_id "FBgnMixed"; gene_symbol "Mixed"; transcript_id "FBtrMixed";',
                        '3R\tFlyBase\tCDS\t350\t390\t.\t-\t0\tgene_id "FBgnMixed"; gene_symbol "Mixed"; transcript_id "FBtrMixed";',
                    ]
                )
                + "\n"
            )
            all_genes = tmp / "all.tsv"
            cds_genes = tmp / "cds.tsv"
            reps = tmp / "reps.tsv"
            cds_gtf = tmp / "cds.gtf"
            exon_gtf = tmp / "exon.gtf"
            rep_cds = tmp / "rep_cds.tsv"
            subprocess.run(
                [
                    "python3",
                    str(ROOT / "workflow/scripts/build_gene_annotation.py"),
                    "--gtf",
                    str(gtf),
                    "--annotation-release",
                    "test_release",
                    "--all-genes",
                    str(all_genes),
                    "--cds-genes",
                    str(cds_genes),
                    "--representative-transcripts",
                    str(reps),
                    "--gene-cds-gtf",
                    str(cds_gtf),
                    "--gene-exon-gtf",
                    str(exon_gtf),
                    "--representative-cds",
                    str(rep_cds),
                ],
                check=True,
            )
            with all_genes.open() as handle:
                gene = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(gene["strand"], ".")
            union_strands = {line.split("\t")[6] for line in cds_gtf.read_text().splitlines()}
            self.assertEqual(union_strands, {"+", "-"})
            with reps.open() as handle:
                rep_rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual([row["transcript_id"] for row in rep_rows], ["FBtrMinus"])


if __name__ == "__main__":
    unittest.main()
