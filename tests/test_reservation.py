from __future__ import annotations

from pathlib import Path
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]


class ReservationTests(unittest.TestCase):
    def test_changed_input_requires_new_analysis_id(self):
        with tempfile.TemporaryDirectory() as tmp_name:
            tmp = Path(tmp_name)
            source = tmp / "decisions.tsv"
            source.write_text("first\n", encoding="utf-8")
            result_dir = tmp / "results/study/a1"
            output = result_dir / "analysis_metadata.json"
            command = [
                "python3", str(ROOT / "workflow/scripts/reserve_analysis.py"),
                "--result-dir", str(result_dir), "--study", "study", "--analysis-id", "a1",
                "--input", str(source), "--output", str(output),
            ]
            subprocess.run(command, check=True)
            subprocess.run(command, check=True)
            source.write_text("changed\n", encoding="utf-8")
            changed = subprocess.run(command, capture_output=True, text=True)
            self.assertNotEqual(changed.returncode, 0)
            self.assertIn("choose a new analysis_id", changed.stderr)


if __name__ == "__main__":
    unittest.main()
