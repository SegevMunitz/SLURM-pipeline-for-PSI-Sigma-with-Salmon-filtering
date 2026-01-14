#!/bin/bash
#SBATCH --job-name=MAKE_GROUPS
#SBATCH --partition=glacier
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:10:00
#SBATCH --output=/sci/labs/zvika.granot/segev.munitz/slurm_logs/%x_%j.out
#SBATCH --error=/sci/labs/zvika.granot/segev.munitz/slurm_logs/%x_%j.err

# -------------------------------------------------------------------
# HOW TO SUBMIT THIS JOB
#
# This is a single (non-array) job.
# It should be run AFTER all STAR alignment array jobs complete.
#
# Submit:
#   sbatch sbatch_03_make_groups.slurm
#
# Or chained after STAR:
#   jid=$(sbatch --array=1-${N}%8 sbatch_02_star_array.slurm | awk '{print $4}')
#   sbatch --dependency=afterok:${jid} sbatch_03_make_groups.slurm
#
# What it does:
#   - Runs step3_make_group_files.py
#   - Generates PSI-Sigma group files from STAR BAM outputs
# -------------------------------------------------------------------

set -euo pipefail

# Ensure we run from the submission directory
cd "$SLURM_SUBMIT_DIR"

python3 step3_make_group_files.py
