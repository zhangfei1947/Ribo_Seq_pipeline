from __future__ import annotations

import os
from pathlib import Path
import subprocess
import tempfile
import textwrap
import unittest


ROOT = Path(__file__).resolve().parents[1]


class ProjectIsolationTests(unittest.TestCase):
    def test_workstation_wrapper_uses_project_working_directory(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_name:
            tmp = Path(tmp_name)
            project_root = tmp / "project" / "analysis"
            config_dir = tmp / "configuration"
            config_dir.mkdir()
            config = config_dir / "config.yaml"
            config.write_text(
                textwrap.dedent(
                    f"""
                    directories:
                      project_root: {project_root}
                      generated: {project_root}/resources/generated
                    """
                ).lstrip(),
                encoding="utf-8",
            )

            fake_bin = tmp / "bin"
            fake_bin.mkdir()
            fake_snakemake = fake_bin / "snakemake"
            fake_snakemake.write_text(
                "#!/usr/bin/env bash\nprintf 'cwd=%s\\n' \"$PWD\"\nprintf 'args=%s\\n' \"$*\"\n",
                encoding="utf-8",
            )
            fake_snakemake.chmod(0o755)

            environment = os.environ.copy()
            environment["PATH"] = f"{fake_bin}:{environment['PATH']}"
            completed = subprocess.run(
                [
                    str(ROOT / "bin" / "run_workstation"),
                    "-n",
                    "qc",
                    "--configfile",
                    str(config.relative_to(tmp)),
                ],
                cwd=tmp,
                env=environment,
                check=True,
                capture_output=True,
                text=True,
            )

            self.assertIn(f"cwd={project_root}", completed.stdout)
            self.assertIn(f"--configfile {config.resolve()}", completed.stdout)
            self.assertTrue((project_root / ".snakemake" / "tmp").is_dir())
            self.assertNotEqual(
                (ROOT / ".snakemake" / "tmp").resolve(),
                (project_root / ".snakemake" / "tmp").resolve(),
            )

    def test_snakefile_has_no_pipeline_global_generated_outputs(self) -> None:
        snakefile = (ROOT / "workflow" / "Snakefile").read_text(encoding="utf-8")
        self.assertNotIn('"resources/generated/', snakefile)
        self.assertIn('GENERATED = config["directories"].get(', snakefile)


if __name__ == "__main__":
    unittest.main()
