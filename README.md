# RNA-seq preprocessing pipeline for Psi-Sigma (SLURM-based)

This repository contains a scalable, generic RNA-seq preprocessing pipeline
designed to generate **STAR alignment outputs suitable for alternative splicing analysis with PSI-Sigma**.

The pipeline is optimized for HPC environments using SLURM arrays and supports
large numbers of samples.

> **Important note**
> If the only goal is **alternative splicing analysis with PSI-Sigma**, you only need:
>
> * **Step 1 – Preparation** (one-time, optional if a STAR index already exists)
> * **Step 2 – STAR alignment** (required)
> * **Step 5 – PSI-Sigma analysis**
>
> Gene-level counting with featureCounts is **optional** and **not required** for PSI-Sigma.

---

## Pipeline overview

The pipeline is composed of the following execution stages:

1. **Preparation** (one-time job; optional if STAR index exists)
2. **STAR alignment** (SLURM array, one sample per task)
3. **FASTQ preprocessing** (optional)
4. **Gene counting with featureCounts** (optional, one-time aggregation job)
5. **PSI-Sigma analysis** (optional downstream splicing analysis)

Each stage is submitted as a separate SLURM job with explicit dependencies.

---

## Directory structure (key files)

* `paths.py`
  Central configuration file (paths, samples, parameters)

* `step1_prep_run.py`
  **Step 1** – Preparation: creates `samples.txt` and builds the STAR index

* `step2_star_align_core.py`
  **Step 2** – STAR alignment worker, designed to run as a SLURM array task

* `step3_fastq_preprocessing.py`
  **Optional Step 3** FASTQ QC and trimming (FastQC, MultiQC, fastp)

* `step4_featurecounts.py`
  **Optional Step 4** – gene-level counting using featureCounts

* `step5_run_psi_sigma.py`
  **Step 5** – PSI-Sigma runner (uses STAR outputs to perform splicing analysis)

* `helper_step0_build_star_index.py`
  STAR genome index builder (used by step1_prep_run.py)

* `sbatch_01_prep.slurm`
  SLURM script for **Step 1 - preparation**

* `sbatch_02_star_array.slurm`
  SLURM array script for **Step 2 – STAR alignment**

* `sbatch_04_featurecounts.slurm`
  SLURM script for **Optional Step 4 – gene-level counting**

* `sbatch_05_run_psi_sigma.slurm`
  SLURM script for **Step 5 – PSI-Sigma analysis**

---

## Step-by-step execution

### Step 1: Preparation (sbatch_01_prep.slurm) — *optional*

Runs `step1_prep_run.py`.

What it does:

* Creates the output directory (`OUTDIR`)
* Writes `OUTDIR/samples.txt` (one sample ID per line)
* Builds the STAR genome index (if not already present)

This step runs **once per genome / annotation / read-length combination**.
If you already have a STAR index and a `samples.txt` file, this step can be skipped.

---

### Step 2: STAR alignment (sbatch_02_star_array.slurm) — **required**

Runs `step2_star_align_core.py` as a SLURM array job.

What it does:

* Each array task processes exactly one sample
* Sample selection is based on `SLURM_ARRAY_TASK_ID`
* FASTQ file paths are read **directly from paths.SAMPLES**
* STAR alignment is performed
* BAM files are coordinate-sorted and indexed
* STAR splice junction files (`SJ.out.tab`) are generated

Output (per sample):

* `OUTDIR/bam/<sample>/<sample>.Aligned.sortedByCoord.out.bam`
* `OUTDIR/bam/<sample>/<sample>.Aligned.sortedByCoord.out.bam.bai`
* `OUTDIR/bam/<sample>/<sample>.SJ.out.tab`

These files are **the required inputs for PSI-Sigma**.

---

### Step 3: FASTQ preprocessing - *optional*

Implemented in `step3_fastq_preprocessing.py`.

Includes:

* FastQC + MultiQC on raw FASTQs
* Optional trimming with fastp
* FastQC + MultiQC on trimmed FASTQs

This step is optional and not required for PSI-Sigma.

---

### Step 4: Gene counting (sbatch_04_featurecounts.slurm) - *optional*

Runs `step4_featurecounts.py` as a single job.

What it does:

* Collects all STAR BAM files
* Runs featureCounts
* Produces a clean gene-level count matrix

Outputs:

* `OUTDIR/counts/featureCounts.gene_counts.txt`
* `OUTDIR/counts/gene_counts.matrix.tsv`

> **Note**
> PSI-Sigma does **not** use gene-level count matrices.
> This step is provided only for users who also want differential gene expression
> or joint expression–splicing analyses.

---

### Step 5: PSI-Sigma analysis (sbatch_05_run_psi_sigma.slurm)

Runs `step5_run_psi_sigma.py` as a single SLURM job.

What it does:

* Uses STAR splice junction files (`SJ.out.tab`) and BAM files
* Applies PSI-Sigma to detect differential alternative splicing events
* Compares predefined sample groups (defined in `paths.py`)

Inputs:

* `OUTDIR/bam/<sample>/<sample>.Aligned.sortedByCoord.out.bam`
* `OUTDIR/bam/<sample>/<sample>.SJ.out.tab`
* Group definition files listed in `paths.py` (e.g. `GROUP_HEALTHY`, `GROUP_SICK_*`)

Outputs:

* PSI-Sigma result tables (event-level splicing statistics)
* Diagnostic and log files produced by PSI-Sigma

---

## How to run this pipeline

This section describes exactly what a user must configure and run in order
to execute the pipeline successfully on a SLURM cluster.

---

## 1. Prerequisites

The following tools must be available on the system (via PATH or modules):

* Python 3.9+
* STAR
* samtools
* SLURM workload manager

Optional (only if running optional steps):

* featureCounts (Subread)
* fastqc
* multiqc
* fastp

---

## 2. Configure paths.py — mandatory

`paths.py` is the single source of configuration for the pipeline.

You must edit the following variables.

### Output

* `OUTDIR`: Root output directory for the entire run

### Reference files

* `GENOME_FASTA`:   Absolute path to the reference genome FASTA file
* `ANNOTATION_GTF`: Absolute path to the gene annotation GTF file
* `STAR_INDEX_DIR`: Directory where the STAR index will be created or reused
* `FASTQ_DIR`: Directory containing the input FASTQ files.
  FASTQ files are discovered automatically from this directory and used to
  construct the `SAMPLES` list. There is no need to manually edit `SAMPLES`.

* `SAMPLES`: A list of dictionaries, one per sample, **generated automatically**
  from the contents of `FASTQ_DIR`.

Each sample dictionary has the following structure:

* `id`: sample name (derived from the FASTQ filename)
* `r1`: absolute path to the R1 FASTQ file
* `r2`: absolute path to the R2 FASTQ file (or omitted for single-end data)

FASTQ filenames must follow the pattern `<sample_id>_R1.fastq.gz` / `<sample_id>_R2.fastq.gz`.
If your naming scheme differs, edit `R1_SUFFIX` and `R2_SUFFIX` in `paths.py`.

Only samples with complete R1/R2 pairs are included.


## 3. Edit SLURM scripts

All SLURM scripts must be edited before running.

### sbatch_01_prep.slurm (Step 0)

* Set partition, memory, and time limits
* Ensure it runs:

  ```bash
  python step1_prep_run.py
  ```

### sbatch_02_star_array.slurm (Step 1)

* Set partition, CPUs, memory, and time
* Set `OUTDIR` and `STAR_INDEX_DIR`
* The array range (`1-N`) is **NOT hardcoded**
* The array range is provided at submission time

### sbatch_04_featurecounts.slurm (Step 2, optional)

* Only needed if gene-level counts are desired

---
## 4. Submit jobs with dependencies

### Minimal PSI-Sigma workflow

```bash
prep_id=$(sbatch sbatch_01_prep.slurm | awk '{print $4}')

N=$(wc -l < "$OUTDIR/samples.txt")
star_id=$(sbatch --array=1-${N}%6 --dependency=afterok:$prep_id sbatch_02_star_array.slurm | awk '{print $4}')

sbatch --dependency=afterok:$star_id sbatch_05_run_psi_sigma.slurm
```

### Optional full workflow (STAR + featureCounts + PSI-Sigma)

```bash
prep_id=$(sbatch sbatch_01_prep.slurm | awk '{print $4}')

N=$(wc -l < "$OUTDIR/samples.txt")
star_id=$(sbatch --array=1-${N}%6 --dependency=afterok:$prep_id sbatch_02_star_array.slurm | awk '{print $4}')

fc_id=$(sbatch --dependency=afterok:$star_id sbatch_04_featurecounts.slurm | awk '{print $4}')

sbatch --dependency=afterok:$star_id sbatch_05_run_psi_sigma.slurm
```

---

## 5. PSI-Sigma outputs

After successful completion of **Step 5**, PSI-Sigma produces:

* Differential splicing event tables
* Event-level PSI statistics
* Summary and log files

Exact output files depend on the PSI-Sigma configuration and group definitions.

---

## 6. Common failure points

* Incorrect paths in `paths.py`
* Sample IDs derived from FASTQ filenames not matching expectations
* STAR index built with the wrong `READ_LENGTH`
* Incorrect SLURM array size
* Missing tools on PATH
* Incorrect group definition files for PSI-Sigma
