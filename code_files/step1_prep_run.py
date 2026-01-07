
"""
step0_prep_run.py

Preparation step for a scalable SLURM-based RNA-seq alignment workflow.

This script performs the one-time setup work needed before launching a STAR
alignment SLURM array. It does NOT run alignment itself.

What it does
------------
1) Creates the pipeline output directory (paths.OUTDIR).
2) Writes a samples list file:
       <OUTDIR>/samples.txt
   containing one sample ID per line, taken from paths.SAMPLES[*]["id"].
   This file is consumed by the STAR SLURM-array worker
   (e.g., step1_star_align_core.py) which uses SLURM_ARRAY_TASK_ID to select
   which sample to process.

3) Builds a STAR genome index (if it does not already exist) by calling
   build_star_index(...) with paths-provided configuration:
     - paths.GENOME_FASTA
     - paths.ANNOTATION_GTF
     - paths.STAR_INDEX_DIR
     - paths.THREADS
     - paths.READ_LENGTH

Expected inputs (from paths.py)
-------------------------------
Required variables in paths.py:
  OUTDIR: Path
      Output root for the run (will be created if missing).

  SAMPLES: list[dict]
      Each dict must include at least:
        {"id": <sample_id>, "r1": <path>, "r2": <path or None>}
      Only the "id" field is used by this script; r1/r2 are used later.

  GENOME_FASTA: Path
      Reference genome FASTA used for STAR index generation.

  ANNOTATION_GTF: Path
      Gene annotation GTF used for STAR index generation.

  STAR_INDEX_DIR: Path
      Directory where the STAR index will be created (or reused).

  THREADS: int
      Threads for STAR genomeGenerate.

  READ_LENGTH: int
      Read length used to compute sjdbOverhang = READ_LENGTH - 1 for index generation.

Outputs
-------
Creates/updates:
  <OUTDIR>/samples.txt
      One sample ID per line.

Creates (if not already present):
  <STAR_INDEX_DIR>/
      STAR index files (sentinel typically <STAR_INDEX_DIR>/Genome).

How this fits in the SLURM workflow
-----------------------------------
Typical run order:
  1) Submit this script as a single SLURM job (prep).
  2) Submit the STAR alignment SLURM array, pointing at:
        --samples <OUTDIR>/samples.txt
  3) After the array completes, submit a final aggregation job (e.g., featureCounts).

Example SLURM usage
-------------------
In a prep SLURM script:
  python step_prep_run.py

Then launch the alignment array:
  python step1_star_align_core.py \
    --samples <OUTDIR>/samples.txt \
    --fastq-dir <FASTQ_DIR> \
    --genome-dir <STAR_INDEX_DIR> \
    --out-dir <OUTDIR>

Notes
-----
- This script assumes the STAR index build step is safe to run multiple times;
  build_star_index(...) should skip if the index already exists.
- This script does not validate FASTQ paths; that validation happens in later steps.
- Ensure the sample IDs in samples.txt match the FASTQ naming convention expected
  by your alignment worker (or pass --r1-suffix/--r2-suffix accordingly).
"""
#!/usr/bin/env python3

from __future__ import annotations
from pathlib import Path

import paths
from utils import ensure_dir, write_text
from helper_step0_build_star_index import build_star_index

def main() -> None:
    out = paths.OUTDIR.resolve()
    ensure_dir(out)

    # samples.txt for the SLURM array worker (one sample ID per line)
    samples_txt = out / "samples.txt"
    write_text(samples_txt, "\n".join([str(s["id"]) for s in paths.SAMPLES]) + "\n")
    print(f"Wrote: {samples_txt}")

    # Build index (one-time)
    build_star_index(
        genome_fasta=paths.GENOME_FASTA,
        annotation_gtf=paths.ANNOTATION_GTF,
        star_index_dir=paths.STAR_INDEX_DIR,
        threads=paths.THREADS,
        read_length=paths.READ_LENGTH,
    )

    print("Prep done.")
    print("Next: submit STAR array using step1_star_align_core.py")

if __name__ == "__main__":
    main()
