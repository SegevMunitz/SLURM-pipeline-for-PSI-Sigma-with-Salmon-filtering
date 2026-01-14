#!/usr/bin/env python3
"""
paths.py

Central configuration and sample discovery for the RNA-seq SLURM pipeline.

This module acts as the single source of truth for paths, resources, and run
settings used across all pipeline steps (prep, STAR alignment, featureCounts,
and PSI-Sigma). Scripts import values from here to avoid duplicating hardcoded
paths in multiple places.

Key responsibilities:
  1) Define reference resources (GENOME_FASTA, ANNOTATION_GTF, STAR_INDEX_DIR).
  2) Define run outputs (OUTDIR) and raw data location (FASTQ_DIR).
  3) Discover samples automatically by scanning FASTQ_DIR and building SAMPLES,
     a list of dicts containing {"id", "r1", "r2"} for downstream steps.
  4) Define optional tool configuration (Perl base, PSI-Sigma script path),
     and group files/comparisons used by PSI-Sigma.

Most users should only edit the variables under the "EDIT THESE PATHS" section.
The remainder of the file is derived configuration used by the pipeline.
"""

from __future__ import annotations
from pathlib import Path

from dataclasses import dataclass

# ----------------- Helpers -----------------

# Local Perl tree (only needed if PSI-Sigma depends on user-local Perl modules like PDL, Statistics::R, etc.)
PERLBASE = Path("/sci/labs/zvika.granot/segev.munitz/softwares/perl5")
# PSI-Sigma executable (dummyai.pl)
PSI_SIGMA = Path("/sci/labs/zvika.granot/segev.munitz/softwares/PSI-Sigma-2.1/dummyai.pl")
# Salmon executable
SALMON_BIN = Path("/sci/labs/zvika.granot/segev.munitz/softwares/Salmon/salmon-latest_linux_x86_64/bin/salmon")

# Probably do not need to edit roots
REFERENCES_ROOT = Path("/sci/labs/zvika.granot/segev.munitz/references")
STAR_INDEX_ROOT  = Path("/sci/labs/zvika.granot/segev.munitz/star_indexes")
SALMON_INDEX_ROOT = Path("/sci/labs/zvika.granot/segev.munitz/salmon_indexes")

R1_SUFFIX = "_1.fastq.gz"
R2_SUFFIX = "_2.fastq.gz"

@dataclass(frozen=True)
class Comparison:
    """
    Describe a single PSI-Sigma comparison to run.

    Each comparison has a human-readable name and two group list files.
    groupa and groupb should be text files containing one absolute BAM path
    per line (created after STAR alignment).
    """
    name: str
    groupa: Path
    groupb: Path

def build_samples_from_fastq_dir(fastq_dir: Path,
                                 r1_suffix: str = R1_SUFFIX,
                                 r2_suffix: str = R2_SUFFIX,
                                 require_paired: bool = True) -> list[dict]:
    """
    Auto-generate SAMPLES by scanning a FASTQ directory for R1/R2 pairs.

    Finds files matching '*<r1_suffix>' and uses the filename prefix as sample ID.
    If require_paired=True, only complete R1/R2 pairs are included; otherwise
    single-end samples are allowed (r2 omitted when absent).
    """
    fastq_dir = Path(fastq_dir)
    if not fastq_dir.is_dir():
        raise SystemExit(f"ERROR: FASTQ_DIR does not exist or is not a directory: {fastq_dir}")

    samples: list[dict] = []

    # Find all R1 files and infer sample IDs
    r1_files = sorted(fastq_dir.glob(f"*{r1_suffix}"))
    if not r1_files:
        raise SystemExit(
            f"ERROR: No R1 files found in {fastq_dir} matching pattern '*{r1_suffix}'. "
            f"Check FASTQ_DIR and suffix."
        )

    for r1 in r1_files:
        sample_id = r1.name[: -len(r1_suffix)]
        r2 = fastq_dir / f"{sample_id}{r2_suffix}"

        if require_paired:
            if not r2.exists():
                # skip incomplete pairs
                continue
            samples.append({"id": sample_id, "r1": r1, "r2": r2})
        else:
            entry = {"id": sample_id, "r1": r1}
            if r2.exists():
                entry["r2"] = r2
            samples.append(entry)

    if require_paired and not samples:
        raise SystemExit(
            f"ERROR: Found R1 files but no complete R1/R2 pairs using "
            f"suffixes r1='{r1_suffix}', r2='{r2_suffix}' in {fastq_dir}."
        )

    # Stable ordering
    samples.sort(key=lambda d: d["id"])
    return samples

#######################################
# --------- EDIT THESE PATHS ----------
#######################################

# Unique name for this pipeline run
NAME = "run_1"

# Pipeline output root
OUTDIR = Path("/sci/labs/zvika.granot/segev.munitz/psi_sigma_outputs") / NAME

# fill FASTQ input dir
FASTQ_DIR = OUTDIR / "gz_files"

# STAR index output dir, Choose Human/Mouse and version as needed
KIND = "Mouse"              # Human or Mouse
VERSION = "M38"             # GENCODE version, e.g. v45 for human, vM31 for mouse
READ_LENGTH = 100           # needed for STAR sjdbOverhang = READ_LENGTH-1

# -------------- DO NOT EDIT -----------------
GENOME_PRE_FIX = "GRCh38" if KIND == "Human" else "GRCm39"      # Reference kind and version
# --------------------------------------------

# Group files: absolute BAM paths, one per line
GROUP_HEALTHY = OUTDIR / "groups" / "H.bams.txt"
GROUP_SICK_1  = OUTDIR / "groups" / "N1.bams.txt"
GROUP_SICK_2  = OUTDIR / "groups" / "N2.bams.txt"
# Add more groups as needed

# Comparisons to run: (name, group_a, group_b)
COMPARISONS: list[Comparison] = [
    Comparison("H_vs_N1", GROUP_HEALTHY, GROUP_SICK_1),
    Comparison("H_vs_N2", GROUP_HEALTHY, GROUP_SICK_2),
    Comparison("N1_vs_N2", GROUP_SICK_1, GROUP_SICK_2),
    # add more as needed
]

# ---------------------------
# GENCODE URLs (optional download)
# ---------------------------

# If you want auto-download to populate REF_BUNDLE_DIR when files are missing:
# NOTE: Human vs Mouse use different FTP paths. Set this correctly for Mouse.
if KIND == "Human":
    # Human release_45 for v45
    GENCODE_BASE_URL = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45"
else:
    # Mouse example (YOU MUST ADJUST to match your VERSION/release on disk)
    # e.g. ".../Gencode_mouse/release_M31" if VERSION="vM31"
    GENCODE_BASE_URL = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M38"

#######################################
# ------ END OF EDITING SECTION -------
#######################################

# Derived reference locations

# Salmon transcripts URL
GENOME_URL = f"{GENCODE_BASE_URL}/{GENOME_PRE_FIX}.primary_assembly.genome.fa.gz"
ANNOTATION_URL = f"{GENCODE_BASE_URL}/gencode.{VERSION}.annotation.gtf.gz"
TRANSCRIPTS_URL = f"{GENCODE_BASE_URL}/gencode.{VERSION}.transcripts.fa.gz"

# Folder name pattern:
#   <root>/<KIND>/Genecode_<VERSION>_rl<READ_LENGTH>/
REF_BUNDLE_DIR = REFERENCES_ROOT / KIND / f"Genecode_{VERSION}_rl{READ_LENGTH}"

# STAR index directory (built files)
STAR_INDEX_DIR = STAR_INDEX_ROOT / KIND / f"Genecode_{VERSION}_rl{READ_LENGTH}"

# Salmon index directory (built files)
SALMON_INDEX_DIR = SALMON_INDEX_ROOT / KIND / f"Genecode_{VERSION}_rl{READ_LENGTH}"

# Reference input files (stored under REFERENCES_ROOT, not inside indexes)
GENOME_FASTA = REF_BUNDLE_DIR / f"{GENOME_PRE_FIX}.primary_assembly.genome.fa"
ANNOTATION_GTF = REF_BUNDLE_DIR / f"gencode.{VERSION}.annotation.gtf"
SALMON_TRANSCRIPTS_FASTA = REF_BUNDLE_DIR / f"gencode.{VERSION}.transcripts.fa"

# FASTQ files are discovered automatically from FASTQ_DIR using suffixes above.
SAMPLES = build_samples_from_fastq_dir(FASTQ_DIR, require_paired=True)

#######################################
# -------=-- OTHER SETTINGS -----------
#######################################

# Compute
THREADS = 16

# Optional steps
DO_QC = False       # FastQC + MultiQC
DO_TRIM = False    # fastp

# featureCounts strandedness: 0=unstranded, 1=stranded, 2=reverse
FEATURECOUNTS_STRANDED = 0

# Extra STAR args if you want
STAR_EXTRA_ARGS = ["--twopassMode", "Basic"] # add more as needed

# Salmon filtering options
SALMON_EXTRA_ARGS = ""
SALMON_LIBTYPE = "A"
SALMON_TPM_THRESHOLD = 5.0   # mean TPM cutoff
SALMON_FILTER_MODE = "either"  # "either" or "both"

PSISIGMA_ABSPSI_MIN: float = 20.0
PSISIGMA_P_MAX: float = 0.05
PSISIGMA_FDR_MAX: float = 0.05

#!/usr/bin/env python3