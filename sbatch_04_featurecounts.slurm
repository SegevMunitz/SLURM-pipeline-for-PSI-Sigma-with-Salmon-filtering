#!/bin/bash
#SBATCH --job-name=FEATURE_COUNTS
#SBATCH --partition=glacier
#SBATCH --cpus-per-task=12
#SBATCH --mem=60G
#SBATCH --time=24:00:00
#SBATCH --output=/sci/labs/zvika.granot/segev.munitz/slurm_logs/%x_%j.out
#SBATCH --error=/sci/labs/zvika.granot/segev.munitz/slurm_logs/%x_%j.err

# -------------------------------------------------------------------
# HOW TO SUBMIT THIS JOB
#
# This is a single (non-array) job.
# It should be run AFTER STAR alignment completes.
#
# Submit:
#   sbatch sbatch_04_featurecounts.slurm
#
# Optional chaining after STAR:
#   jid=$(sbatch --array=1-${N}%8 sbatch_02_star_array.slurm | awk '{print $4}')
#   sbatch --dependency=afterok:${jid} sbatch_04_featurecounts.slurm
#
# What it does:
#   - Runs step4_featurecounts.py
#   - Generates gene-level read counts using featureCounts
#
# Options:
#   --paired-end : use paired-end counting mode
# -------------------------------------------------------------------

set -euo pipefail

# Ensure we run from the submission directory
cd "$SLURM_SUBMIT_DIR"

python3 step4_featurecounts.py \
    --paired-end
