# RNA-seq SLURM Pipeline: STAR Alignment & PSI-Sigma Splicing

This repository provides a scalable, modular pipeline for processing RNA-seq data on HPC clusters using the SLURM scheduler. It is specifically optimized to generate **STAR alignment outputs and differential alternative splicing results with PSI-Sigma**.

## Table of Contents
1. Pipeline Overview
2. Prerequisites
3. Configuration & Setup
4. Execution Steps
5. Automated Job Submission (Dependencies)
6. Output Directory Structure

---

## Pipeline Overview
The pipeline is divided into five distinct stages, each paired with a Python "core" script and a corresponding .slurm submission script:

**Step	        Script	                    Description**
1. **Prep**	    `step1_prep_run.py`	        Builds the STAR index and initializes the directory structure.
2. **Align**	`step2_star_align_core.py`	Aligns FASTQ files to the genome using SLURM job arrays.
3. **Group**	`step3_make_group_files.py`	Organizes BAM files into experimental groups (e.g., Control vs. Treatment).
4. **Count**	`step4_featurecounts.py`	(Optional) Generates a gene expression matrix.
5. **Splicing**	`step5_run_psi_sigma.py`  	Executes PSI-Sigma to find differential splicing events.

---

## Prerequisites & Environment
Ensure the following are available in your `$PATH`:
* Python 3.8+
* STAR (v2.7.x recommended)
* Samtools (for BAM indexing)
* Subread (for featureCounts)
* Perl (with PDL and Statistics::R modules for PSI-Sigma)
* PSI-Sigma (the dummyai.pl script and its dependencies)
---

## Configuration & Setup
The pipeline uses a "Single Source of Truth" configuration file: `paths.py`.

1. Edit `paths.py`:
Open `paths.py` and update the following critical variables:

* `GENOME_FASTA / ANNOTATION_GTF`: Reference genome and gene model.
* `OUTDIR`: The root directory for all pipeline outputs.
* `FASTQ_DIR`: Location of your .fastq.gz files.
* `COMPARISONS`: Define your experimental pairs (e.g., Healthy vs Sick).

2. Update SLURM Scripts
In each .slurm file, update the following under the # EDIT section:

* `#SBATCH --output and --error: Ensure these point to a valid log directory`.
* `#SBATCH --partition: Match your cluster's partition name (e.g., glacier)`.

---

## Execution Steps
**Step 1: Initialization & Indexing**
Builds the STAR index if missing and creates `samples.txt` file used for the array job.

```bash
sbatch sbatch_01_prep.slurm
```

---

**Step 2: STAR Alignment (Array)**
This is the most resource-intensive step. It launches a parallel task for every sample discovered in your FASTQ directory.

```bash
# Calculate the number of samples
N=$(wc -l < /path/to/your/output/samples.txt)

# Submit the array
sbatch --array=1-$N sbatch_02_star_array.slurm
```

> Note: You must set the array range to match the number of lines in samples.txt.

---

**Step 3: Experimental Grouping**
Before running splicing analysis, the pipeline needs to know which samples belong to which group. You must provide a TSV file mapping sample IDs to group names.

```bash
sbatch sbatch_03_make_groups.slurm
```

---

**Step 4: Gene Counting (Optional)**
Runs `featureCounts` across all generated BAM files to create a raw count matrix and a gene-level count matrix for optional differential expression analysis (not used by PSI-Sigma).
```bash
sbatch sbatch_04_featurecounts.slurm
```

---

**Step 5: PSI-Sigma Splicing**
Runs the differential splicing analysis. This script handles the complex Perl environment and temporary file management required by PSI-Sigma.
```bash
sbatch sbatch_05_run_psi_sigma.slurm
```

---

## Automated Job Submission (Dependencies)
To run the entire pipeline without manual intervention between steps, use SLURM's dependency system:

```bash
# 1. Start Prep
PREP_ID=$(sbatch sbatch_01_prep.slurm | awk '{print $4}')

# 2. Get sample count and start STAR (waits for Prep)
# Note: Manually verify N if running for the first time
N=$(ls -1 /path/to/fastq/*_R1.fastq.gz | wc -l)
STAR_ID=$(sbatch --array=1-${N}%6 --dependency=afterok:$PREP_ID sbatch_02_star_array.slurm | awk '{print $4}')

# 3. Splicing Analysis (waits for STAR)
sbatch --dependency=afterok:$STAR_ID sbatch_05_run_psi_sigma.slurm
```

---

## Output Directory Structure
The pipeline maintains a strict, predictable hierarchy:

```
star_indexes/                           # STAR Genome Index
    ├── Human                           # Human Genome Index
        └── Genecode_v45_rl100          # Format for human genomes
        └── ...
    └── Mouse                           # Mouse Genome Index
        └── Genecode_m38_rl100          # Format for mouse genomes
        └── ...
OUTDIR/
├── samples.txt                         # List of sample IDs
├── groups.tsv                          # User-provided mapping
├── bam/                                # Per-sample alignment results
│   └── Sample_01/
│       ├── Sample_01.Aligned.sortedByCoord.out.bam
│       └── Sample_01.Aligned.sortedByCoord.out.bam.bai
├── groups/                             # BAM path lists for PSI-Sigma
│   ├── Healthy.bams.txt
│   └── Sick.bams.txt
├── counts/                             # featureCounts matrices
│   └── gene_counts.matrix.tsv
├── logs/                               # SLURM stdout/stderr
└── Comparison_Name/                    # PSI-Sigma splicing results
    ├── PSIsigma_r10_ir3.sorted.txt     # Filtered Results (p<0.01)
    ├── PSIsigma_r10_ir3.txt            # Unfiltered Results
    ├── *.SJ.out.tab                    # Junction Read File
    ├── *.IR.out.tab                    # Intergenic Read File  
    ├── *.db                            # PSI-Sigma Database
    ├── *.bam                           # BAM File
    └── *.gtf                           # GTF File
```

---

## Troubleshooting
* STAR Index Failures: Ensure you have enough disk space and memory (at least 32GB for human/mouse genomes). Check `helper_step0_build_star_index.py` logs.
* Sample Discovery: If samples aren't being found, check `R1_SUFFIX` and `R2_SUFFIX` in `paths.py`. They must exactly match your filename endings.
* PSI-Sigma Perl Errors: If PSI-Sigma fails due to missing Perl modules (like PDL), ensure `PERLBASE` in `paths.py` points to a Perl installation where those modules are installed.
* SLURM Array Limits: Some clusters limit the number of simultaneous tasks. The `%6` in `sbatch_02_star_array.slurm` limits the array to 6 concurrent jobs; adjust this based on your cluster's policy.
