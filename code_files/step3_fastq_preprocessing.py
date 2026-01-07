
"""
step2_fastq_preprocessing.py

Utilities for FASTQ quality control (QC) and optional read trimming prior to
downstream analysis (e.g., STAR alignment, Salmon, kallisto).

This module provides a small, explicit preprocessing layer that can be inserted
at the beginning of an RNA-seq pipeline. It supports:
  - FastQC + MultiQC on raw FASTQs
  - Optional adapter/quality trimming with fastp
  - FastQC + MultiQC on trimmed FASTQs
  - Returning the FASTQ paths that should be used downstream

The goal is to make QC and trimming decisions explicit, reproducible, and easy
to enable or disable.

Overview
--------
Given a list of samples (each with an ID and FASTQ paths), the typical flow is:

  FASTQ
    → FastQC + MultiQC (raw)
    → [optional] fastp trimming
    → [optional] FastQC + MultiQC (trimmed)
    → FASTQs ready for alignment / quantification

No alignment or quantification is performed in this module.

Functions
---------
1) fastqc_multiqc(...)
   Runs FastQC on a list of FASTQ files and aggregates the results with MultiQC.

2) fastp_trim(...)
   Runs fastp on one sample (single-end or paired-end), producing trimmed FASTQs
   and fastp HTML/JSON reports.

3) process_fastqs(...)
   High-level orchestrator that:
     - validates FASTQ paths
     - runs QC on raw FASTQs (optional)
     - trims reads with fastp (optional)
     - runs QC on trimmed FASTQs (optional)
     - returns the FASTQ paths to use downstream

Key design decisions
--------------------
- QC and trimming are optional and controlled by flags (`do_qc`, `do_trim`).
- Trimming is NOT automatic; it must be explicitly enabled.
- QC on trimmed reads is only run if trimming was performed.
- All external tools (fastqc, multiqc, fastp) are checked via `which_or_die`
  before execution.
- The module returns FASTQ paths rather than writing global state, allowing
  clean integration into downstream pipeline stages.

Input data model
----------------
Samples are provided as a list of dictionaries with the following keys:

  {
    "id": "sample_name",
    "r1": "/path/to/sample_R1.fastq.gz",
    "r2": "/path/to/sample_R2.fastq.gz",   # optional (omit or set to None for SE)
  }

Single-end samples are supported by omitting the "r2" field.

Output layout
-------------
All outputs are written under `out_root`:

  out_root/
    qc_raw/            FastQC + MultiQC reports for raw FASTQs (if enabled)
    qc_trimmed/        FastQC + MultiQC reports for trimmed FASTQs (if enabled)
    trimmed_fastq/     Trimmed FASTQs and fastp reports (if trimming enabled)

Trimming outputs per sample:
  <sample>.R1.trimmed.fastq.gz
  <sample>.R2.trimmed.fastq.gz   (paired-end only)
  <sample>.fastp.html
  <sample>.fastp.json

Return value
------------
process_fastqs(...) returns a dictionary:

  {
    sample_id: (Path_to_R1, Path_to_R2_or_None),
    ...
  }

If trimming is disabled, the returned paths point to the original FASTQs.
If trimming is enabled, the returned paths point to the trimmed FASTQs.

When to use this module
----------------------
- At the start of a new RNA-seq analysis
- When working with new sequencing data or a new sequencing facility
- When you want explicit QC reports for documentation or publication
- When you want trimming to be reproducible and optional

When NOT to use trimming
------------------------
- If FASTQs were already trimmed by the sequencing provider
- If FastQC on raw data shows clean adapter content and high base quality
- If you need strict comparability with a previous analysis that did not trim

Examples
--------
Basic usage with QC only:
  trimmed = process_fastqs(
      samples=samples,
      out_root=Path("qc"),
      threads=8,
      do_qc=True,
      do_trim=False,
  )

QC + trimming:
  trimmed = process_fastqs(
      samples=samples,
      out_root=Path("qc"),
      threads=8,
      do_qc=True,
      do_trim=True,
  )

Downstream integration:
  for sid, (r1, r2) in trimmed.items():
      run_star(sample_id=sid, r1=r1, r2=r2)

Notes
-----
- This module does not modify FASTQs in place.
- STAR and other modern aligners tolerate minor low-quality bases, so trimming
  should be driven by QC results, not habit.
- All heavy lifting is delegated to external, well-established tools.
"""

from __future__ import annotations
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from utils import ensure_dir, run, which_or_die

def fastqc_multiqc(fastqs: List[Path], out_dir: Path, threads: int, label: str) -> None:
    which_or_die("fastqc")
    which_or_die("multiqc")
    qc_dir = out_dir / f"qc_{label}"
    ensure_dir(qc_dir)
    run(["fastqc", "-t", str(threads), "-o", str(qc_dir)] + [str(fq) for fq in fastqs])
    run(["multiqc", str(qc_dir), "-o", str(qc_dir)])

def fastp_trim(
    sample_id: str,
    r1: Path,
    r2: Optional[Path],
    out_dir: Path,
    threads: int,
) -> Tuple[Path, Optional[Path]]:
    which_or_die("fastp")
    ensure_dir(out_dir)

    out_r1 = out_dir / f"{sample_id}.R1.trimmed.fastq.gz"
    out_r2 = out_dir / f"{sample_id}.R2.trimmed.fastq.gz" if r2 else None
    html = out_dir / f"{sample_id}.fastp.html"
    js = out_dir / f"{sample_id}.fastp.json"

    cmd = ["fastp", "-w", str(threads), "-i", str(r1), "-o", str(out_r1), "-h", str(html), "-j", str(js)]
    if r2:
        cmd += ["-I", str(r2), "-O", str(out_r2)]
    run(cmd)

    return out_r1, out_r2

def process_fastqs(
    samples: List[Dict],
    out_root: Path,
    threads: int,
    do_qc: bool,
    do_trim: bool,
) -> Dict[str, Tuple[Path, Optional[Path]]]:
    # Validate inputs
    for s in samples:
        if not Path(s["r1"]).exists():
            raise SystemExit(f"ERROR: missing R1 for {s['id']}: {s['r1']}")
        if s.get("r2") and not Path(s["r2"]).exists():
            raise SystemExit(f"ERROR: missing R2 for {s['id']}: {s['r2']}")

    # QC raw
    if do_qc:
        raw_fastqs: List[Path] = []
        for s in samples:
            raw_fastqs.append(Path(s["r1"]))
            if s.get("r2"):
                raw_fastqs.append(Path(s["r2"]))
        fastqc_multiqc(raw_fastqs, out_root, threads, "raw")

    # Trim (or passthrough)
    trimmed: Dict[str, Tuple[Path, Optional[Path]]] = {}
    trim_dir = out_root / "trimmed_fastq"

    for s in samples:
        sid = s["id"]
        r1 = Path(s["r1"])
        r2 = Path(s["r2"]) if s.get("r2") else None

        if do_trim:
            trimmed[sid] = fastp_trim(sid, r1, r2, trim_dir, threads)
        else:
            trimmed[sid] = (r1, r2)

    # QC trimmed
    if do_qc and do_trim:
        t_fastqs: List[Path] = []
        for sid, (r1, r2) in trimmed.items():
            t_fastqs.append(r1)
            if r2:
                t_fastqs.append(r2)
        fastqc_multiqc(t_fastqs, out_root, threads, "trimmed")

    return trimmed
