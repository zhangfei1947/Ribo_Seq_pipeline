#!/usr/bin/env python3
"""Translate a SLURM job ID into Snakemake's generic status vocabulary."""

from __future__ import annotations

import subprocess
import sys


RUNNING = {"PENDING", "CONFIGURING", "RUNNING", "COMPLETING", "SUSPENDED", "REQUEUED"}
SUCCESS = {"COMPLETED"}
FAILED = {
    "BOOT_FAIL",
    "CANCELLED",
    "DEADLINE",
    "FAILED",
    "NODE_FAIL",
    "OUT_OF_MEMORY",
    "PREEMPTED",
    "REVOKED",
    "TIMEOUT",
}


def query(job_id: str) -> str:
    sacct = subprocess.run(
        ["sacct", "-j", job_id, "--format=State", "--noheader", "--parsable2"],
        check=False,
        capture_output=True,
        text=True,
    )
    states = [line.strip().split()[0].split("+")[0] for line in sacct.stdout.splitlines() if line.strip()]
    if not states:
        queue = subprocess.run(
            ["squeue", "-h", "-j", job_id, "-o", "%T"],
            check=False,
            capture_output=True,
            text=True,
        )
        states = [line.strip().split()[0] for line in queue.stdout.splitlines() if line.strip()]
    if any(state in FAILED for state in states):
        return "failed"
    if states and all(state in SUCCESS for state in states):
        return "success"
    return "running"


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("usage: slurm_status.py JOB_ID")
    print(query(sys.argv[1]))

