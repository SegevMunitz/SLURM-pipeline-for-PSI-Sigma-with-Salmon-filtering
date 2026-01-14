# RNA-seq SLURM Pipeline: STAR Alignment, PSI-Sigma Splicing & Salmon Quantification

This repository provides a scalable, modular pipeline for processing RNA-seq data on HPC clusters using the SLURM scheduler. It is specifically optimized to generate **STAR alignment outputs, differential alternative splicing results with PSI-Sigma, and transcript quantification with Salmon**.

## Table of Contents
1. [Pipeline Overview](#pipeline-overview)
2. [Prerequisites & Environment](#prerequisites--environment)
3. [Configuration & Setup](#configuration--setup)
4. [Execution Steps](#execution-steps)
5. [Automated Job Submission](#automated-job-submission-dependencies)
6. [Output Directory Structure](#output-directory-structure)
7. [Troubleshooting](#troubleshooting)

---

## Pipeline Overview

The pipeline is divided into **seven distinct stages**, each paired with a Python "core" script and a corresponding `.slurm` submission script:

| Step | Script | Description |
|------|--------|-------------|
| **1. Prep** | `step1_prep_run.py` | Builds STAR and Salmon indexes, downloads references if needed, and initializes directory structure |
| **2. Align** | `step2_star_align_core.py` | Aligns FASTQ files to the genome using STAR via SLURM job arrays |
| **3. Group** | `step3_make_group_files.py` | Organizes BAM files into experimental groups based on `groups.tsv` |
| **4. Count** | `step4_featurecounts.py` | Generates gene expression count matrix using featureCounts |
| **5. Splicing** | `step5_run_psi_sigma.py` | Executes PSI-Sigma to find differential splicing events |
| **6. Salmon** | `step6_salmon_quant_core.py` | Quantifies transcript expression using Salmon via SLURM job arrays |
| **7. Filter** | `step7_psi_sigma_filtered_by_salmon.py` | Filters PSI-Sigma results based on Salmon TPM thresholds |

---

## Prerequisites & Environment

Ensure the following tools are available in your `$PATH`:

### Required Software
- **Python 3.8+**
- **STAR** (v2.7.x recommended) - RNA-seq aligner
- **Samtools** - for BAM indexing
- **Subread** (featureCounts) - for gene counting
- **Salmon** - for transcript quantification
- **Perl** (with PDL and Statistics::R modules) - for PSI-Sigma
- **PSI-Sigma** (the `dummyai.pl` script and dependencies)

### Python Packages
```bash
pip install pandas openpyxl --break-system-packages
```

---

## Configuration & Setup

The pipeline uses a **"Single Source of Truth"** configuration file: `paths.py`.

### 1. Edit `paths.py`

Open `paths.py` and update the following critical variables:

#### Core Paths (EDIT THESE)
```python
# Unique name for this pipeline run
NAME = "run_1"

# Pipeline output root
OUTDIR = Path("/path/to/your/output") / NAME

# FASTQ input directory
FASTQ_DIR = OUTDIR / "gz_files"  # or specify a different path

# Genome selection
KIND = "Mouse"              # "Human" or "Mouse"
VERSION = "M38"             # GENCODE version (e.g., "v45" for Human, "M38" for Mouse)
READ_LENGTH = 100           # Your sequencing read length

# Group files for PSI-Sigma comparisons
GROUP_HEALTHY = OUTDIR / "groups" / "H.bams.txt"
GROUP_SICK_1  = OUTDIR / "groups" / "N1.bams.txt"
GROUP_SICK_2  = OUTDIR / "groups" / "N2.bams.txt"

# Define your comparisons
COMPARISONS: list[Comparison] = [
    Comparison("H_vs_N1", GROUP_HEALTHY, GROUP_SICK_1),
    Comparison("H_vs_N2", GROUP_HEALTHY, GROUP_SICK_2),
    Comparison("N1_vs_N2", GROUP_SICK_1, GROUP_SICK_2),
]
```

#### File Naming Convention
The pipeline expects paired-end FASTQ files with these suffixes (configurable):
```python
R1_SUFFIX = "_1.fastq.gz"  # e.g., sample1_1.fastq.gz
R2_SUFFIX = "_2.fastq.gz"  # e.g., sample1_2.fastq.gz
```

#### Software Paths
Update these paths to match your installation:
```python
PERLBASE = Path("/path/to/perl5")  # Local Perl installation with PDL/Statistics::R
PSI_SIGMA = Path("/path/to/PSI-Sigma-2.1/dummyai.pl")
SALMON_BIN = Path("/path/to/salmon/bin/salmon")
```

#### Reference Paths
The pipeline automatically constructs reference paths based on `KIND`, `VERSION`, and `READ_LENGTH`:
```python
REFERENCES_ROOT = Path("/path/to/references")
STAR_INDEX_ROOT = Path("/path/to/star_indexes")
SALMON_INDEX_ROOT = Path("/path/to/salmon_indexes")
```

This creates:
- Reference files: `REFERENCES_ROOT/<KIND>/Genecode_<VERSION>_rl<READ_LENGTH>/`
- STAR index: `STAR_INDEX_ROOT/<KIND>/Genecode_<VERSION>_rl<READ_LENGTH>/`
- Salmon index: `SALMON_INDEX_ROOT/<KIND>/Genecode_<VERSION>_rl<READ_LENGTH>/`

#### Optional Settings
```python
THREADS = 16                          # Default threads for indexing
FEATURECOUNTS_STRANDED = 0            # 0=unstranded, 1=stranded, 2=reverse
STAR_EXTRA_ARGS = ["--twopassMode", "Basic"]  # Additional STAR parameters

# Salmon filtering
SALMON_TPM_THRESHOLD = 5.0            # Minimum mean TPM for filtering
SALMON_FILTER_MODE = "either"         # "either" or "both" groups must pass threshold

# PSI-Sigma filtering
PSISIGMA_ABSPSI_MIN: float = 20.0     # Minimum |ΔPSI| threshold
PSISIGMA_P_MAX: float = 0.05          # Maximum p-value
PSISIGMA_FDR_MAX: float = 0.05        # Maximum FDR
```

### 2. Update SLURM Scripts

In each `.slurm` file in `sbatch_commands/`, update the following:

```bash
#SBATCH --output=/path/to/logs/%x_%A_%a.out
#SBATCH --error=/path/to/logs/%x_%A_%a.err
#SBATCH --partition=your_partition_name
```

Also ensure the scripts point to the correct working directory containing your Python scripts.

### 3. Create `groups.tsv`

Create a file `<OUTDIR>/groups.tsv` mapping samples to groups:

```tsv
sample1	H
sample2	H
sample3	N1
sample4	N1
sample5	N2
sample6	N2
```

Format: `<sample_id><TAB><group_name>` (no header)

**Important:** Sample IDs must match the filenames in your FASTQ directory (without suffixes).

---

## Execution Steps

### Step 1: Initialization & Indexing

Builds STAR and Salmon indexes if missing, downloads reference files if needed, and creates `samples.txt`.

```bash
cd sbatch_commands
sbatch sbatch_01_prep.slurm
```

**What it does:**
- Downloads genome FASTA, GTF, and transcripts if not present
- Builds STAR genome index (if missing)
- Builds Salmon transcript index (if missing)
- Creates `samples.txt` from FASTQ directory scan

---

### Step 2: STAR Alignment (Array Job)

This is the most resource-intensive step. It launches a parallel task for every sample discovered in your FASTQ directory.

```bash
# Count samples
N=$(wc -l < <OUTDIR>/samples.txt)

# Submit the array (limit to 8 concurrent jobs with %8)
sbatch --array=1-$N%8 sbatch_02_star_array.slurm
```

**What it does:**
- Runs STAR alignment for each sample
- Outputs sorted BAM files
- Creates BAM indexes (.bai)
- Runs in two-pass mode for improved junction detection

**Output location:** `<OUTDIR>/star_alignments/bam/<sample>/<sample>.Aligned.sortedByCoord.out.bam`

---

### Step 3: Experimental Grouping

Creates group files containing BAM paths for PSI-Sigma analysis.

```bash
sbatch sbatch_03_make_groups.slurm
```

**Requirements:**
- Must have `<OUTDIR>/groups.tsv` file (see Configuration above)
- All BAM files from Step 2 must exist

**What it does:**
- Reads `groups.tsv` mapping
- Validates that all samples have BAM files in `<OUTDIR>/star_alignments/bam/`
- Creates `<OUTDIR>/groups/<group_name>.bams.txt` files with absolute BAM paths

---

### Step 4: Gene Counting

Runs `featureCounts` to create a gene expression matrix.

```bash
sbatch sbatch_04_featurecounts.slurm
```

**What it does:**
- Auto-discovers BAMs under `<OUTDIR>/bam/*/*.Aligned.sortedByCoord.out.bam` (or uses explicit `--bam` paths)
- Counts reads per gene across all samples
- Generates raw featureCounts output
- Creates cleaned matrix: `<OUTDIR>/counts/gene_counts.matrix.tsv`

**Note:** If BAMs are in `<OUTDIR>/star_alignments/bam/`, you may need to either:
- Create a symlink: `ln -s star_alignments/bam <OUTDIR>/bam`
- Or provide explicit `--bam` arguments to the script

**Note:** This step is independent from PSI-Sigma but useful for differential expression analysis.

---

### Step 5: PSI-Sigma Splicing Analysis

Runs differential splicing analysis for all comparisons defined in `paths.py`.

```bash
sbatch sbatch_05_run_psi_sigma.slurm
```

**What it does:**
- Runs PSI-Sigma for each comparison in `COMPARISONS`
- Detects differential splicing events
- Outputs results per comparison

**Output:** `<OUTDIR>/<comparison_name>/PSIsigma_*.txt`

---

### Step 6: Salmon Quantification (Array Job)

Quantifies transcript-level expression for all samples.

```bash
# Use same N from Step 2
sbatch --array=1-$N%8 sbatch_06_salmon_array.slurm
```

**What it does:**
- Runs Salmon on FASTQ files for each sample
- Produces transcript TPM values
- Handles paired-end and single-end data

**Output:** `<OUTDIR>/salmon/<sample>/quant.sf`

---

### Step 7: Filter PSI-Sigma by Salmon Expression

Filters PSI-Sigma results based on transcript expression levels.

```bash
sbatch sbatch_07_salmon_filter.slurm
```

**What it does:**
- Loads Salmon TPM values for each group
- Adds mean TPM columns to PSI-Sigma results
- Filters by expression threshold and statistical significance
- Exports filtered results as Excel and TSV

**Output:** `<OUTDIR>/<comparison_name>/<comparison>.PSI<threshold>_p<pval>_FDR<fdr>.salmon_tpm<tpm>.<mode>.xlsx`

---

## Automated Job Submission (Dependencies)

To run the entire pipeline without manual intervention, use the provided submission script:

```bash
cd sbatch_commands
chmod +x submit_all.sh
./submit_all.sh
```

Or manually chain jobs with dependencies:

```bash
# 1. Prep
jid_prep=$(sbatch sbatch_01_prep.slurm | awk '{print $4}')

# 2. Count samples
N=$(wc -l < <OUTDIR>/samples.txt)

# 3. STAR alignment (after prep)
jid_star=$(sbatch --dependency=afterok:${jid_prep} \
                  --array=1-${N}%8 \
                  sbatch_02_star_array.slurm | awk '{print $4}')

# 4. Make groups (after STAR)
jid_groups=$(sbatch --dependency=afterok:${jid_star} \
                    sbatch_03_make_groups.slurm | awk '{print $4}')

# 5. featureCounts (after STAR)
jid_fc=$(sbatch --dependency=afterok:${jid_star} \
                sbatch_04_featurecounts.slurm | awk '{print $4}')

# 6. PSI-Sigma (after groups + featureCounts)
jid_psisigma=$(sbatch --dependency=afterok:${jid_groups}:${jid_fc} \
                      sbatch_05_run_psi_sigma.slurm | awk '{print $4}')

# 7. Salmon (after STAR)
jid_salmon=$(sbatch --dependency=afterok:${jid_star} \
                    --array=1-${N}%8 \
                    sbatch_06_salmon_array.slurm | awk '{print $4}')

# 8. Filter (after PSI-Sigma + Salmon)
sbatch --dependency=afterok:${jid_psisigma}:${jid_salmon} \
       sbatch_07_salmon_filter.slurm
```

**Dependency Graph:**
```
Prep (1)
  ↓
STAR Array (2)
  ↓
  ├─→ Groups (3) ──→ PSI-Sigma (5) ──┐
  ├─→ featureCounts (4) ─────────────┤
  └─→ Salmon Array (6) ──────────────┤
                                      ↓
                              Filter Results (7)
```

---

## Output Directory Structure

```
<OUTDIR>/
├── samples.txt                           # List of sample IDs (created by step1, used by step2)
├── groups.tsv                            # User-provided sample-to-group mapping (you create this)
│
├── star_alignments/                      # STAR alignment outputs (created by step2)
│   ├── bam/                             # BAM files (accessed by step3)
│   │   └── <sample>/
│   │       ├── <sample>.Aligned.sortedByCoord.out.bam
│   │       ├── <sample>.Aligned.sortedByCoord.out.bam.bai
│   │       ├── <sample>.Log.final.out
│   │       ├── <sample>.Log.out
│   │       ├── <sample>.Log.progress.out
│   │       └── <sample>.SJ.out.tab     # Splice junctions
│   ├── logs/                            # STAR-specific logs
│   └── tmp/                             # STAR temporary files
│
├── groups/                               # BAM path lists for PSI-Sigma (created by step3)
│   ├── H.bams.txt
│   ├── N1.bams.txt
│   └── N2.bams.txt
│
├── counts/                               # featureCounts outputs (created by step4)
│   ├── featureCounts.gene_counts.txt    # Raw output
│   └── gene_counts.matrix.tsv           # Cleaned matrix
│
├── salmon/                               # Salmon quantification (created by step6)
│   └── <sample>/
│       ├── quant.sf                     # Transcript TPMs
│       ├── quant.genes.sf               # Gene-level aggregation
│       ├── lib_format_counts.json
│       └── aux_info/
│
├── <comparison_name>/                    # PSI-Sigma results (created by step5, filtered by step7)
│   ├── PSIsigma_r10_ir3.txt            # Unfiltered results
│   ├── <comparison>.PSI20_p0.05_FDR0.05.salmon_tpm5.0.either.xlsx  # Filtered by step7
│   └── <comparison>.PSI20_p0.05_FDR0.05.salmon_tpm5.0.either.tsv   # Filtered by step7
│
├── logs/                                 # General pipeline logs
└── tmp/                                  # General temporary files
```

### Reference Directories (External to OUTDIR)

```
<REFERENCES_ROOT>/<KIND>/Genecode_<VERSION>_rl<READ_LENGTH>/
├── GRCh38.primary_assembly.genome.fa     # Genome FASTA (Human)
├── GRCm39.primary_assembly.genome.fa     # Genome FASTA (Mouse)
├── gencode.v45.annotation.gtf            # Annotation GTF (Human)
├── gencode.M38.annotation.gtf            # Annotation GTF (Mouse)
└── gencode.<VERSION>.transcripts.fa      # Transcripts for Salmon

<STAR_INDEX_ROOT>/<KIND>/Genecode_<VERSION>_rl<READ_LENGTH>/
├── Genome
├── SA
├── SAindex
└── [other STAR index files]

<SALMON_INDEX_ROOT>/<KIND>/Genecode_<VERSION>_rl<READ_LENGTH>/
├── complete_ref_lens.bin
├── ctable.bin
├── ctg_offsets.bin
├── duplicate_clusters.tsv
├── info.json
└── [other Salmon index files]
```

---

## Troubleshooting

### STAR Index Failures
- **Issue:** Index building fails or runs out of memory
- **Solution:** 
  - Ensure at least 32GB RAM for human/mouse genomes
  - Check `helper_step0_build_star_index.py` logs
  - Verify genome FASTA and GTF are properly downloaded and uncompressed

### Sample Discovery Issues
- **Issue:** No samples found or wrong samples detected
- **Solution:**
  - Check that `R1_SUFFIX` and `R2_SUFFIX` in `paths.py` exactly match your filenames
  - Example: If files are `sample1_R1.fq.gz`, set `R1_SUFFIX = "_R1.fq.gz"`
  - Verify FASTQ files are in the correct `FASTQ_DIR`
  - Check that R1 and R2 pairs exist for all samples (if `require_paired=True`)

### PSI-Sigma Perl Errors
- **Issue:** Missing Perl modules (PDL, Statistics::R, etc.)
- **Solution:**
  - Ensure `PERLBASE` in `paths.py` points to a valid Perl installation with required modules
  - Check Perl environment: `perl -MPDL -e 'print "OK\n"'`
  - Install missing modules: `cpanm PDL Statistics::R`

### featureCounts Cannot Find BAM Files
- **Issue:** Step 4 reports "No BAMs found"
- **Root cause:** featureCounts auto-discovery looks for `<OUTDIR>/bam/*/*.bam`, but BAMs are in `<OUTDIR>/star_alignments/bam/`
- **Solution 1 (Symlink):**
  ```bash
  ln -s star_alignments/bam <OUTDIR>/bam
  ```
- **Solution 2 (Explicit paths):**
  Modify the sbatch script to pass explicit `--bam` arguments for each sample

### Missing Group Files
- **Issue:** Step 3 fails with "groups.tsv not found"
- **Solution:**
  - Create `<OUTDIR>/groups.tsv` before running Step 3
  - Format: `<sample_id><TAB><group_name>` (no header)
  - Ensure sample IDs match those in `samples.txt`

### Salmon Quantification Fails
- **Issue:** Salmon index not found or quantification errors
- **Solution:**
  - Verify Salmon index was built successfully in Step 1
  - Check `SALMON_BIN` path in `paths.py`
  - Ensure transcripts FASTA was downloaded correctly
  - For array jobs: Only task 1 builds index, others wait (this is intentional)

### Array Job Limits
- **Issue:** Too many jobs pending or cluster policy violations
- **Solution:**
  - The `%8` in array submissions limits concurrent jobs to 8
  - Adjust this based on your cluster policy: `--array=1-$N%<max_concurrent>`
  - Check with `squeue -u $USER` to monitor job status

### File Not Found Errors
- **Issue:** Pipeline can't find expected BAM or output files
- **Solution:**
  - Verify previous steps completed successfully
  - Check SLURM log files in your specified log directory
  - Ensure file paths in `paths.py` are absolute, not relative
  - Check disk space and quotas: `df -h` and `quota -s`

### Memory or Timeout Errors
- **Issue:** Jobs fail with out-of-memory or timeout errors
- **Solution:**
  - Increase memory in SLURM script: `#SBATCH --mem=60G` → `#SBATCH --mem=120G`
  - Increase time limit: `#SBATCH --time=24:00:00` → `#SBATCH --time=48:00:00`
  - Reduce `--outBAMsortingThreadN` in `step2_star_align_core.py` if memory limited

### Salmon Filtering Produces No Results
- **Issue:** Step 7 outputs empty files
- **Solution:**
  - Lower `SALMON_TPM_THRESHOLD` in `paths.py` (e.g., from 5.0 to 1.0)
  - Change `SALMON_FILTER_MODE` from "both" to "either"
  - Verify Salmon quantification produced valid TPM values in `quant.sf` files
  - Check that Reference_Transcript IDs in PSI-Sigma match Salmon transcript IDs

### Permission Denied Errors
- **Issue:** Cannot write to output directories
- **Solution:**
  - Check directory permissions: `ls -ld <OUTDIR>`
  - Ensure you have write access to all paths in `paths.py`
  - Create directories manually if needed: `mkdir -p <OUTDIR>`

---

## Additional Notes

### Customizing STAR Parameters
To modify STAR behavior, edit these variables in `paths.py`:
```python
STAR_EXTRA_ARGS = [
    "--twopassMode", "Basic",
    "--outFilterMultimapNmax", "20",
    "--alignSJoverhangMin", "8",
    # Add more parameters as needed
]
```

### Single-End Data
The pipeline supports single-end reads. Simply ensure:
- `require_paired=False` in `build_samples_from_fastq_dir()` in `paths.py`
- Only R1 files in FASTQ_DIR
- Remove `-p -B -C` flags from featureCounts in Step 4

### Rerunning Failed Steps
If a step fails:
1. Fix the issue (check logs in `<OUTDIR>/logs/`)
2. Resubmit only that step's SLURM script
3. Dependencies will automatically trigger downstream steps

### Contact & Support
For issues related to:
- **PSI-Sigma:** [PSI-Sigma GitHub](https://github.com/wososa/PSI-Sigma)
- **STAR:** [STAR GitHub](https://github.com/alexdobin/STAR)
- **Salmon:** [Salmon documentation](https://salmon.readthedocs.io/)

---

**Pipeline Version:** 2.0  
**Last Updated:** January 2025

