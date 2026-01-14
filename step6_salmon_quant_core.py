#!/usr/bin/env python3
"""
step6_salmon_quant_core.py

Run Salmon quantification on all samples using a SLURM array job.
This script is intended to be run as part of a SLURM array job, where each
task quantifies one sample.
"""

import subprocess
import time
from pathlib import Path
import os
import sys

from paths import (
    SAMPLES,
    OUTDIR,
    SALMON_BIN,
    SALMON_INDEX_DIR,
    SALMON_LIBTYPE,
    SALMON_EXTRA_ARGS,
    SALMON_TRANSCRIPTS_FASTA,
    TRANSCRIPTS_URL,
    GENOME_FASTA,
)

# Fix helper script path if needed
SALMON_HELPER = Path(__file__).parent / "helper_step1_build_salmon_index.py"

def salmon_index_exists(index_dir: Path) -> bool:
    """
    Check if a valid Salmon index exists.
    
    A minimal reliable check: both info.json and ctable.bin must be present.
    """
    return (index_dir / "info.json").exists() and (index_dir / "ctable.bin").exists()


def ensure_salmon_index():
    """
    Build Salmon index if missing.
    
    In SLURM array jobs, only task 1 builds the index.
    Other tasks wait until the index appears.
    
    This prevents race conditions where multiple tasks try to build
    the same index simultaneously.
    """
    if salmon_index_exists(SALMON_INDEX_DIR):
        return

    task_id = int(os.environ.get("SLURM_ARRAY_TASK_ID", "1"))

    if task_id == 1:
        print("Salmon index not found - building now")

        cmd = [
            sys.executable,
            str(SALMON_HELPER),
            "--transcripts-fasta", str(SALMON_TRANSCRIPTS_FASTA),
            "--index-dir", str(SALMON_INDEX_DIR),
            "--threads", "16",
            "--download-dir", str(SALMON_TRANSCRIPTS_FASTA.parent),
            "--transcripts-url", str(TRANSCRIPTS_URL),
            "--genome-fasta", str(GENOME_FASTA)
        ]
        subprocess.run(cmd, check=True)

    else:
        print("Waiting for Salmon index to be built by task 1...")
        start = time.time()
        timeout_s = 6 * 60 * 60  # 6 hours

        while not salmon_index_exists(SALMON_INDEX_DIR):
            if time.time() - start > timeout_s:
                raise SystemExit(f"ERROR: Timed out waiting for Salmon index at {SALMON_INDEX_DIR}")
            time.sleep(30)


def run_salmon(sample_id, fastq1, fastq2=None):
    """
    Run Salmon quantification for a single sample.
    
    Args:
        sample_id: Sample identifier (used for output directory)
        fastq1: Path to R1 FASTQ file (or single-end FASTQ)
        fastq2: Path to R2 FASTQ file (None for single-end)
    """
    outdir = OUTDIR / "salmon" / sample_id
    outdir.mkdir(parents=True, exist_ok=True)

    cmd = [
        SALMON_BIN, "quant",
        "-i", str(SALMON_INDEX_DIR),
        "-l", SALMON_LIBTYPE,
        "-o", str(outdir),
        "--gcBias",
        "--validateMappings",
    ]

    if SALMON_EXTRA_ARGS:
        cmd.extend(SALMON_EXTRA_ARGS.split())

    if fastq2:
        # Paired-end
        cmd.extend(["-1", fastq1, "-2", fastq2])
    else:
        # Single-end
        cmd.extend(["-r", fastq1])

    subprocess.run(cmd, check=True)


if __name__ == "__main__":

    # Ensure index exists before any quantification
    ensure_salmon_index()

    task_id = int(os.environ["SLURM_ARRAY_TASK_ID"]) - 1
    if task_id >= len(SAMPLES):
        raise SystemExit(f"ERROR: SLURM_ARRAY_TASK_ID={task_id+1} exceeds number of samples ({len(SAMPLES)})")
    sample = SAMPLES[task_id]

    sample_id = sample["id"]
    fq1 = str(sample["r1"])
    fq2 = str(sample.get("r2")) if "r2" in sample else None

    run_salmon(sample_id, fq1, fq2)

    print(f"\nSalmon quantification complete for {sample_id}")
