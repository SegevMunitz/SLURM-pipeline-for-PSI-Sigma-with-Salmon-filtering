#!/usr/bin/env python3

"""
step1_star_align.py

Run STAR alignment for paired-end FASTQ samples using an existing STAR genome index.
Designed to mirror a typical SLURM-array bash workflow.

This script does NOT build the index. Use step0_get_star_index.py to build the index first.

Core behavior (per selected sample)
-----------------------------------
1) Selects a sample name from a samples file (one sample per line).
   - Uses --task-id if provided (1-based).
   - Otherwise uses SLURM_ARRAY_TASK_ID (1-based).
2) Constructs paired-end FASTQ paths:
       R1 = <fastq-dir>/<sample><r1-suffix>
       R2 = <fastq-dir>/<sample><r2-suffix>
   Defaults:
       r1-suffix = "_1.fastq.gz"
       r2-suffix = "_2.fastq.gz"
3) Creates output directories:
       <out-dir>/logs/
       <out-dir>/bam/<sample>/
       <out-dir>/tmp/
4) Creates a unique temp directory per run and passes it to STAR:
       <out-dir>/tmp/<sample>.tmp_<jobid>_<taskid>
   If it already exists, it is removed and recreated (like `rm -rf` in bash).
5) Runs STAR with parameters matching the bash defaults:
       --runThreadN <threads>
       --outSAMtype BAM SortedByCoordinate
       --outFilterIntronMotifs RemoveNoncanonical
       --outBAMsortingThreadN 4
       --genomeDir <index-dir>
       --twopassMode Basic
       --readFilesCommand zcat        (optional; can be disabled)
       --readFilesIn <R1> <R2>
       --outFileNamePrefix <out-dir>/bam/<sample>/<sample>.
       --outTmpDir <tmp-dir>
6) Verifies the expected BAM exists:
       <out-dir>/bam/<sample>/<sample>.Aligned.sortedByCoord.out.bam
7) Indexes the BAM:
       samtools index <bam>

Thread selection
----------------
- If --threads is provided, it is used.
- Else SLURM_CPUS_PER_TASK is used if set.

CLI arguments
-------------
Required:
  --genome-dir PATH
  --fastq-dir PATH
  --out-dir PATH
  --samples PATH

Optional:
  --task-id INT          1-based line number in samples file
  --threads INT
  --r1-suffix STR
  --r2-suffix STR
  --tmp-base PATH

STAR parameter overrides (optional):
  --twopassMode STR
  --outFilterIntronMotifs STR
  --outBAMsortingThreadN INT
  --readFilesCommand STR   Set to empty string to disable (for uncompressed FASTQ)

Examples
--------
SLURM array usage (task id inferred from environment):
  python star_align.py \\
    --genome-dir /indexes/mm10_rl100 \\
    --fastq-dir /data/fastq \\
    --out-dir /results/star_out \\
    --samples /results/star_out/samples.txt

Local single task:
  python star_align.py \\
    --genome-dir /indexes/mm10_rl100 \\
    --fastq-dir /data/fastq \\
    --out-dir /results/star_out \\
    --samples /results/star_out/samples.txt \\
    --task-id 3 \\
    --threads 12

Uncompressed FASTQ:
  python star_align.py \\
    --genome-dir /indexes/mm10_rl100 \\
    --fastq-dir /data/fastq \\
    --out-dir /results/star_out \\
    --samples /results/star_out/samples.txt \\
    --readFilesCommand ""
"""

from __future__ import annotations

from pathlib import Path
import argparse
import os
import socket

from paths import SAMPLES
from utils import ensure_dir, run, which_or_die

def get_sample_info(sample_id: str):
    for s in SAMPLES:
        if s["id"] == sample_id:
            return s
    raise SystemExit(
        f"ERROR: sample '{sample_id}' not found in paths.SAMPLES "
        f"(available: {[x['id'] for x in SAMPLES]})"
    )

def read_sample_from_list(samples_file: Path, task_id_1based: int) -> str:
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
    which_or_die("samtools")
    run(["samtools", "index", str(bam)])


def build_tmp_dir(tmp_base: Path, sample: str, job_id: str, task_id: str) -> Path:
    tmp = tmp_base / f"{sample}.tmp_{job_id}_{task_id}"
    if tmp.exists():
        # be safe like your bash
        import shutil
        shutil.rmtree(tmp)
    ensure_dir(tmp)
    return tmp


def main() -> None:
    p = argparse.ArgumentParser(description="Run STAR alignment per-sample.")

    # Core IO
    p.add_argument("--genome-dir", type=Path, required=True, help="STAR index directory.")
    p.add_argument("--out-dir", type=Path, required=True)
    p.add_argument("--samples", type=Path, required=True, help="Text file with one sample per line.")

    # Sample/task selection
    p.add_argument("--task-id", type=int, default=None, help="1-based line number in --samples. If omitted, uses SLURM_ARRAY_TASK_ID.")
    p.add_argument("--r1-suffix", type=str, default="_1.fastq.gz")
    p.add_argument("--r2-suffix", type=str, default="_2.fastq.gz")

    # Threads
    p.add_argument("--threads", type=int, default=None, help="If none given, uses SLURM_CPUS_PER_TASK.")

    # STAR params matching bash defaults, chnage if wanted other then default
    p.add_argument("--twopassMode", type=str, default="Basic")
    p.add_argument("--outFilterIntronMotifs", type=str, default="RemoveNoncanonical")
    p.add_argument("--outBAMsortingThreadN", type=int, default=4)
    p.add_argument("--readFilesCommand", type=str, default="zcat", help="Set to '' to disable (for uncompressed FASTQ).")

    # Temp + prefix layout
    p.add_argument("--tmp-base", type=Path, default=None, help="Base tmp directory. Default: <out-dir>/tmp")
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

    out_dir: Path = args.out_dir
    logs_dir = out_dir / "logs"
    bam_root = out_dir / "bam"
    tmp_base = args.tmp_base or (out_dir / "tmp")

    ensure_dir(logs_dir)
    ensure_dir(bam_root)
    ensure_dir(tmp_base)

    # Logging banner (imitate bash)
    job_id = os.environ.get("SLURM_JOB_ID", "NOJOB")
    array_id = os.environ.get("SLURM_ARRAY_TASK_ID", str(task_id))
    host = socket.gethostname()
    print(f"=== job={job_id} task={array_id} host={host} ===")    
    print(f"STAR={star_cmd}")
    print(f"IDX={args.genome_dir}")
    print(f"OUT={out_dir}")
    print(f"SAMPLES={args.samples}")
    print(f"Threads={threads}")

    sample = read_sample_from_list(args.samples, task_id)

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
        "--genomeDir", str(args.genome_dir),
        "--twopassMode", args.twopassMode,
        "--readFilesIn", str(r1),
        "--outFileNamePrefix", str(prefix),
        "--outTmpDir", str(tmp),
    ]

    if r2 is not None:
        cmd.insert(cmd.index("--outFileNamePrefix"), str(r2))  # place r2 right after r1

    if args.readFilesCommand.strip():
        cmd.insert(cmd.index("--readFilesIn"), "--readFilesCommand")
        cmd.insert(cmd.index("--readFilesIn"), args.readFilesCommand)

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
