from __future__ import annotations

import csv
from pathlib import Path
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
FIELDS = ["read_length", "offset", "n_representative_cds", "frame0", "frame1", "frame2", "frame0_fraction"]


class PoolOffsetTests(unittest.TestCase):
    def test_frame_evidence_is_summed_before_selection(self):
        with tempfile.TemporaryDirectory() as tmp_name:
            tmp = Path(tmp_name)
            inputs = []
            for index in range(2):
                path = tmp / f"sample{index}.tsv"
                with path.open("w", newline="", encoding="utf-8") as handle:
                    writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t")
                    writer.writeheader()
                    writer.writerow({"read_length": 28, "offset": 12, "n_representative_cds": 60, "frame0": 48, "frame1": 6, "frame2": 6, "frame0_fraction": 0.8})
                    writer.writerow({"read_length": 28, "offset": 13, "n_representative_cds": 60, "frame0": 20, "frame1": 20, "frame2": 20, "frame0_fraction": 1/3})
                inputs.append(path)
            output = tmp / "pooled.tsv"
            command = ["python3", str(ROOT / "workflow/scripts/pool_offset_recommendations.py")]
            for path in inputs:
                command.extend(["--candidate", str(path)])
            command.extend(["--min-reads-per-length", "100", "--output", str(output)])
            subprocess.run(command, check=True)
            with output.open(encoding="utf-8") as handle:
                row = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(row["recommended_offset"], "12")
            self.assertEqual(row["status"], "recommended")


if __name__ == "__main__":
    unittest.main()
