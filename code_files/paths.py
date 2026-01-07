from __future__ import annotations
from pathlib import Path

from pytoolconfig import dataclass

# ----------------- Sample discovery helpers -----------------

R1_SUFFIX = "_R1.fastq.gz"
R2_SUFFIX = "_R2.fastq.gz"

@dataclass(frozen=True)
class Comparison:
    name: str
    groupa: Path
    groupb: Path

def build_samples_from_fastq_dir(fastq_dir: Path, r1_suffix: str = R1_SUFFIX, r2_suffix: str = R2_SUFFIX, require_paired: bool = True) -> list[dict]:
    """
    Build SAMPLES = [{"id": ..., "r1": Path(...), "r2": Path(...)}] automatically
    by scanning FASTQ_DIR for files matching <id><r1_suffix>.

    - If require_paired=True: only includes samples that have BOTH R1 and R2.
    - If require_paired=False: will include single-end samples too (r2 omitted).
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


# --------- EDIT THESE PATHS ----------
## Choose Human/Mouse
GENOME_FASTA = Path("/sci/labs/zvika.granot/segev_munitz/genome_references/          .fa") # reference genome fasta
ANNOTATION_GTF = Path("/sci/labs/zvika.granot/segev_munitz/genome_references/       .gtf") # annotation gtf
STAR_INDEX_DIR = Path("/sci/labs/zvika.granot/segev_munitz/                             ") # STAR index output dir
OUTDIR = Path("/sci/labs/zvika.granot/segev_munitz/                                     ") # Pipeline output root
FASTQ_DIR = Path("/sci/labs/zvika.granot/segev.munitz/neutrophil_RNAseq/fastq_files/    ") # FASTQ input dir

# FASTQ files are discovered automatically from FASTQ_DIR using suffixes above. No need to edit below.
SAMPLES = build_samples_from_fastq_dir(FASTQ_DIR, require_paired=True)

# Local Perl tree (only needed if PSI-Sigma depends on user-local Perl modules like PDL, Statistics::R, etc.)
PERLBASE = Path("/sci/labs/zvika.granot/segev.munitz/softwares/perl5")

# PSI-Sigma executable (dummyai.pl)
PSI_SIGMA = Path("/sci/labs/zvika.granot/segev.munitz/softwares/PSI-Sigma-2.1/dummyai.pl")

# Group files: absolute BAM paths, one per line

GROUP_HEALTHY = OUTDIR / "groups" / "Healthy.bams.txt"
GROUP_SICK_1  = OUTDIR / "groups" / "Sick_1.bams.txt"
GROUP_SICK_2  = OUTDIR / "groups" / "Sick_2.bams.txt"

# Add more groups as needed

# Compute
THREADS = 16
READ_LENGTH = 100  # needed for STAR sjdbOverhang = READ_LENGTH-1

# Optional steps
DO_QC = False       # FastQC + MultiQC
DO_TRIM = False    # fastp

# featureCounts strandedness: 0=unstranded, 1=stranded, 2=reverse
FEATURECOUNTS_STRANDED = 0

# Extra STAR args if you want
STAR_EXTRA_ARGS = ["--twopassMode", "Basic"]

# Comparisons to run: (name, group_a, group_b)
COMPARISONS: list[Comparison] = [
    Comparison("H_vs_S1", GROUP_HEALTHY, GROUP_SICK_1),
    Comparison("H_vs_S2", GROUP_HEALTHY, GROUP_SICK_2),
    Comparison("S1_vs_S2", GROUP_SICK_1, GROUP_SICK_2),
    # add more as needed
]
