#!/usr/bin/env python3
"""
step1_prep_run.py

Preparation and initialization step for a SLURM-based RNA-seq alignment pipeline.

This script performs all one-time setup tasks required before launching a
STAR alignment SLURM array. It is responsible for preparing the output
directory structure, defining the set of samples to be processed, and
ensuring that the required STAR genome index exists.

Specifically, this script:
  1) Creates the pipeline output root directory defined in paths.OUTDIR.
  2) Writes a canonical samples list file:
         <OUTDIR>/samples.txt
     containing one sample ID per line. This file defines the iteration
     space for the downstream SLURM array via SLURM_ARRAY_TASK_ID.
  3) Invokes STAR genome index generation using configuration provided in
     paths.py. If an index already exists, the build step is skipped.

All configuration is centralized in paths.py, allowing downstream pipeline
steps to rely on a consistent definition of sample identifiers, reference
files, and output locations.

Expected configuration (paths.py)
---------------------------------
Required variables:

  OUTDIR : Path
      Root output directory for the pipeline run.

  SAMPLES : list[dict]
      List of sample metadata dictionaries. Each entry must contain an "id"
      field. FASTQ paths are not validated here and are used in later steps.

  GENOME_FASTA : Path
      Reference genome FASTA for STAR index generation.

  ANNOTATION_GTF : Path
      Gene annotation GTF for STAR index generation.

  STAR_INDEX_DIR : Path
      Directory where the STAR genome index will be created or reused.

  THREADS : int
      Number of threads to use for STAR genome index generation.

  READ_LENGTH : int
      Sequencing read length used to compute sjdbOverhang (READ_LENGTH - 1).

Generated outputs
-----------------
  <OUTDIR>/
      Pipeline output root directory.

  <OUTDIR>/samples.txt
      One sample ID per line, consumed by the STAR SLURM array worker.

  <STAR_INDEX_DIR>/
      STAR genome index files, created if not already present.

Workflow context
----------------
Typical execution order:

  1) Run this script once to prepare the analysis.
  2) Submit the STAR alignment SLURM array, using samples.txt to distribute work.
  3) Run downstream aggregation and analysis steps after alignment completes.

This separation ensures that reference preparation and sample bookkeeping
are performed once, while alignment and downstream processing scale cleanly
across samples using SLURM arrays.
"""

from __future__ import annotations

from paths import (OUTDIR, SAMPLES, STAR_INDEX_DIR, GENOME_FASTA,
                   ANNOTATION_GTF, THREADS, READ_LENGTH, REF_BUNDLE_DIR,
                   GENOME_URL, ANNOTATION_URL, SALMON_INDEX_DIR, SALMON_BIN,
                   TRANSCRIPTS_URL, SALMON_TRANSCRIPTS_FASTA)
from utils import ensure_dir, write_text
from helper_step0_build_star_index import build_star_index
from helper_step1_build_salmon_index import build_salmon_index

def main() -> None:
    """
    Prepare the pipeline output directory and reference resources.
    
    Creates the run output directory, writes samples.txt for downstream
    SLURM array jobs, and ensures a STAR genome index exists by invoking
    the reference build step if needed.
    """
    out = OUTDIR.resolve()
    ensure_dir(out)

    # samples.txt for the SLURM array worker (one sample ID per line)
    samples_txt = out / "samples.txt"
    write_text(samples_txt, "\n".join([str(s["id"]) for s in SAMPLES]) + "\n")
    print(f"Wrote: {samples_txt}")

    # Build STAR index (one-time)
    build_star_index(
        genome_fasta=GENOME_FASTA,
        annotation_gtf=ANNOTATION_GTF,
        star_index_dir=STAR_INDEX_DIR,
        threads=THREADS,
        read_length=READ_LENGTH,
        download_dir=REF_BUNDLE_DIR,
        genome_url=GENOME_URL,
        annotation_url=ANNOTATION_URL
    )

    # Build Salmon index (one-time)
    build_salmon_index(
        transcripts_fasta=SALMON_TRANSCRIPTS_FASTA,
        salmon_index_dir=SALMON_INDEX_DIR,
        threads=THREADS,
        download_dir=REF_BUNDLE_DIR,
        transcripts_url=TRANSCRIPTS_URL
    )

    print("Prep done.")
    print("Next: submit STAR array using step1_star_align_core.py")

if __name__ == "__main__":
    main()
