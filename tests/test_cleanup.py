from __future__ import annotations

import csv
from pathlib import Path
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]


class CleanupTests(unittest.TestCase):
    def test_dry_run_and_confirmation_guard(self):
        with tempfile.TemporaryDirectory() as tmp_name:
            root = Path(tmp_name)
            target = root / "work/alignment/sample.bam"
            target.parent.mkdir(parents=True)
            target.write_bytes(b"bam")
            manifest = root / "manifest.tsv"
            with manifest.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=["relative_path", "bytes", "category", "analysis_id"],
                    delimiter="\t",
                )
                writer.writeheader()
                writer.writerow(
                    {
                        "relative_path": "work/alignment/sample.bam",
                        "bytes": 3,
                        "category": "alignment",
                        "analysis_id": "a1",
                    }
                )
            command = [
                "python3",
                str(ROOT / "workflow/scripts/cleanup_intermediates.py"),
                "--project-root",
                str(root),
                "--manifest",
                str(manifest),
                "--analysis-id",
                "a1",
            ]
            subprocess.run(command, check=True, capture_output=True, text=True)
            self.assertTrue(target.exists())
            failed = subprocess.run(command + ["--execute"], capture_output=True, text=True)
            self.assertNotEqual(failed.returncode, 0)
            subprocess.run(command + ["--execute", "--confirm-analysis-id", "a1"], check=True)
            self.assertFalse(target.exists())


if __name__ == "__main__":
    unittest.main()
