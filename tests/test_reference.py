from __future__ import annotations

import json
from pathlib import Path
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]


class ReferenceValidationTests(unittest.TestCase):
    def make_star(self, path: Path, overhang: int) -> None:
        path.mkdir()
        for name in ("Genome", "SA", "SAindex"):
            (path / name).write_bytes(b"index")
        (path / "chrName.txt").write_text("X\n", encoding="utf-8")
        (path / "genomeParameters.txt").write_text(
            f"sjdbOverhang {overhang}\ngenomeSAindexNbases 12\n", encoding="utf-8"
        )

    def test_compatible_reference(self):
        with tempfile.TemporaryDirectory() as tmp_name:
            tmp = Path(tmp_name)
            ribo, rna = tmp / "ribo", tmp / "rna"
            self.make_star(ribo, 34)
            self.make_star(rna, 149)
            prefix = tmp / "contaminants"
            for suffix in ("1", "2", "3", "4", "rev.1", "rev.2"):
                Path(f"{prefix}.{suffix}.bt2").write_bytes(b"index")
            output = tmp / "report.json"
            subprocess.run(
                [
                    "python3", str(ROOT / "workflow/scripts/validate_reference.py"),
                    "--gtf", str(ROOT / "tests/fixtures/mini.gtf"),
                    "--ribo-star", str(ribo), "--rna-star", str(rna),
                    "--contaminant-prefix", str(prefix), "--ribo-overhang", "34",
                    "--rna-overhang", "149", "--sa-index-bases", "12",
                    "--release", "test", "--output", str(output),
                ],
                check=True,
            )
            self.assertEqual(json.loads(output.read_text())["status"], "valid")


if __name__ == "__main__":
    unittest.main()

