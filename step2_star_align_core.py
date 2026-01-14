#!/usr/bin/env python3
"""
step2_star_align_core.py

Per-sample STAR alignment worker for SLURM arrays.

This script is designed to be executed as part of a SLURM array job, where each
array task aligns exactly one sample selected from a samples list file
(e.g., <OUTDIR>/samples.txt). The selected sample ID is determined by either:
  - --task-id (1-based), or
  - SLURM_ARRAY_TASK_ID (1-based) if --task-id is not provided.

Sample definitions are loaded from paths.SAMPLES, where each entry must include:
  - "id": unique sample identifier (must match a line in samples.txt)
  - "r1": path to R1 FASTQ
  - "r2": optional path to R2 FASTQ (may be None for single-end)

Outputs are written under the provided --out-dir:
  - <OUTDIR>/bam/<sample>/<sample>.Aligned.sortedByCoord.out.bam
  - <OUTDIR>/bam/<sample>/<sample>.Aligned.sortedByCoord.out.bam.bai
  - <OUTDIR>/logs/ (directory created; STAR logs are written in sample folder via prefix)
  - <OUTDIR>/tmp/<sample>.tmp_<jobid>_<taskid> (unique temp per task)

STAR is executed with parameters chosen to match typical pipeline defaults:
  - Sorted BAM output (SortedByCoordinate)
  - Two-pass mode (Basic)
  - Noncanonical intron motif filtering
  - Optional readFilesCommand (zcat by default; auto-disabled if R1 is not .gz)

This script expects STAR and samtools to be available on PATH (module load or conda
activation should happen in the sbatch script). It performs basic validation
(FASTQ existence; BAM creation) and fails fast if required inputs/outputs are missing.
"""
from __future__ import annotations

from pathlib import Path
import argparse
import os
import socket

from paths import (SAMPLES, STAR_INDEX_DIR, OUTDIR)
from utils import ensure_dir, run, which_or_die

def get_sample_info(sample_id: str):
    """
    Look up a sample record in paths.SAMPLES by its sample ID.

    Returns the full sample metadata dict (including r1/r2 paths) used by the
    alignment worker. Exits with a clear error if the sample ID is not found.
    """
    for s in SAMPLES:
        if s["id"] == sample_id:
            return s
    raise SystemExit(
        f"ERROR: sample '{sample_id}' not found in paths.SAMPLES "
        f"(available: {[x['id'] for x in SAMPLES]})"
    )

def read_sample_from_list(samples_file: Path, task_id_1based: int) -> str:
    """
    Read the sample ID corresponding to a 1-based task index from a samples file.

    Mimics the SLURM-array pattern `sed -n "${SLURM_ARRAY_TASK_ID}p"`, selecting
    exactly one line and returning it as the sample ID. Validates index bounds
    and rejects empty/blank lines to prevent silent misalignment.
    """
    if task_id_1based < 1:
        raise SystemExit(f"ERROR: task id must be >= 1 (got {task_id_1based})")
    lines = samples_file.read_text().splitlines()
    # mimic `sed -n "${SLURM_ARRAY_TASK_ID}p"`
    if task_id_1based > len(lines):
        raise SystemExit(f"ERROR: task {task_id_1based} out of range (file has {len(lines)} lines): {samples_file}")
    sample = lines[task_id_1based - 1].strip()
    sample = sample.replace("\r", "").strip()
    if not sample:
        raise SystemExit(f"ERROR: empty sample name at line {task_id_1based} in {samples_file}")
    return sample


def samtools_index(bam: Path) -> None:
    """
    Create a BAM index (.bai) for the given BAM file using samtools.

    Verifies samtools is available on PATH and runs `samtools index`.
    This is typically required for downstream tools (e.g., PSI-Sigma) that
    expect indexed BAMs.
    """
    which_or_die("samtools")
    run(["samtools", "index", str(bam)])


def build_tmp_dir(tmp_base: Path, sample: str, job_id: str, task_id: str) -> Path:
    """
    Create a unique per-task STAR temporary directory under the output root.

    The temp path encodes sample, SLURM job id, and task id to avoid collisions.
    If the directory already exists (e.g., rerun), it is removed and recreated
    to mimic a clean `rm -rf` behavior before STAR runs.
    """
    tmp = tmp_base / f"{sample}.tmp_{job_id}_{task_id}"
    if tmp.exists():
        # be safe like your bash
        import shutil
        shutil.rmtree(tmp)
    ensure_dir(tmp)
    return tmp


def main() -> None:
    """
    Execute STAR alignment for a single sample (one SLURM array task).

    Parses CLI arguments, resolves the task/sample ID, validates FASTQ inputs,
    constructs STAR command-line parameters, runs STAR, verifies the BAM output,
    and indexes the BAM using samtools for downstream analysis.
    """
    p = argparse.ArgumentParser(description="Run STAR alignment per-sample.")

    # Core IO
    genome_dir: Path = STAR_INDEX_DIR
    out_dir: Path = OUTDIR / "star_alignments"
    samples: Path = OUTDIR / "samples.txt"

    # Sample/task selection
    p.add_argument("--task-id", type=int, default=None, help="1-based line number in --samples. If omitted, uses SLURM_ARRAY_TASK_ID.")
    p.add_argument("--r1-suffix", type=str, default="_1.fastq.gz")
    p.add_argument("--r2-suffix", type=str, default="_2.fastq.gz")

    # Threads
    p.add_argument("--threads", type=int, default=None, help="If none given, uses SLURM_CPUS_PER_TASK.")

    # STAR params matching bash defaults, change if wanted other then default
    p.add_argument("--twopassMode", type=str, default="Basic")
    p.add_argument("--outFilterIntronMotifs", type=str, default="RemoveNoncanonical")
    p.add_argument("--outBAMsortingThreadN", type=int, default=4)
    p.add_argument("--readFilesCommand", type=str, default="zcat", help="Set to '' to disable (for uncompressed FASTQ).")

    args = p.parse_args()

    # Resolve STAR executable
    which_or_die("STAR")
    star_cmd = "STAR"

    # Resolve task id
    task_id = args.task_id
    if task_id is None:
        env_tid = os.environ.get("SLURM_ARRAY_TASK_ID")
        if not env_tid:
            raise SystemExit("ERROR: --task-id not provided and SLURM_ARRAY_TASK_ID not set")
        task_id = int(env_tid)

    # Resolve threads
    threads = args.threads
    if threads is None:
        threads = int(os.environ.get("SLURM_CPUS_PER_TASK", "4"))

    logs_dir = out_dir / "logs"
    bam_root = out_dir / "bam"
    tmp_base = out_dir / "tmp"

    ensure_dir(logs_dir)
    ensure_dir(bam_root)
    ensure_dir(tmp_base)

    # Logging banner (imitate bash)
    job_id = os.environ.get("SLURM_JOB_ID", "NOJOB")
    array_id = os.environ.get("SLURM_ARRAY_TASK_ID", str(task_id))
    host = socket.gethostname()
    print(f"=== job={job_id} task={array_id} host={host} ===")    
    print(f"STAR={star_cmd}")
    print(f"IDX={genome_dir}")
    print(f"OUT={out_dir}")
    print(f"SAMPLES={samples}")
    print(f"Threads={threads}")

    sample = read_sample_from_list(samples, task_id)

    rec = get_sample_info(sample)
    r1 = Path(rec["r1"])
    r2_val = rec.get("r2", None)
    r2 = Path(r2_val) if r2_val else None

    if not r1.is_file():
        raise SystemExit(f"ERROR missing R1 for sample={sample}: {r1}")
    if r2 is not None and not r2.is_file():
        raise SystemExit(f"ERROR missing R2 for sample={sample}: {r2}")
    if args.readFilesCommand == "zcat" and not str(r1).endswith(".gz"):
        args.readFilesCommand = ""   # disable for plain fastq

    sample_out = bam_root / sample
    ensure_dir(sample_out)

    prefix = sample_out / f"{sample}."
    tmp = build_tmp_dir(tmp_base, sample, job_id, str(task_id))

    # Run STAR  
    cmd = [
        star_cmd,
        "--runThreadN", str(threads),
        "--outSAMtype", "BAM", "SortedByCoordinate",
        "--outFilterIntronMotifs", args.outFilterIntronMotifs,
        "--outBAMsortingThreadN", str(args.outBAMsortingThreadN),
        "--genomeDir", str(genome_dir),
        "--twopassMode", args.twopassMode,
        "--outFileNamePrefix", str(prefix),
        "--outTmpDir", str(tmp),
    ]
    
    read_files = [str(r1)]

    if r2 is not None:
        read_files.append(str(r2))
    
    cmd += ["--readFilesIn", *read_files]

    if args.readFilesCommand.strip():
        i = cmd.index("--readFilesIn")
        cmd[i:i] = ["--readFilesCommand", args.readFilesCommand]
        
    print(f"--- START STAR sample={sample} ---")
    run(cmd)
    print(f"--- END STAR sample={sample} ---")

    bam = sample_out / f"{sample}.Aligned.sortedByCoord.out.bam"
    if not bam.is_file():
        raise SystemExit(f"ERROR: BAM not created: {bam}")

    print(f"--- START samtools index sample={sample} ---")
    samtools_index(bam)
    print(f"--- END samtools index sample={sample} ---")

    print(f"=== DONE sample={sample} ===")


if __name__ == "__main__":
    main()
