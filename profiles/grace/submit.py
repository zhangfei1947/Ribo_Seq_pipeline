#!/usr/bin/env python3
"""Submit with --parsable so cluster-generic records a numeric SLURM job ID."""

from __future__ import annotations

from pathlib import Path
import subprocess
import sys


def main() -> int:
    Path("logs/slurm").mkdir(parents=True, exist_ok=True)
    completed = subprocess.run(
        ["sbatch", "--parsable", *sys.argv[1:]],
        check=True,
        capture_output=True,
        text=True,
    )
    # Federated SLURM can return JOBID;CLUSTER. Status commands need JOBID only.
    job_id = completed.stdout.strip().split(";", 1)[0]
    if not job_id.isdigit():
        raise SystemExit(f"Unexpected sbatch --parsable output: {completed.stdout!r}")
    print(job_id)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

