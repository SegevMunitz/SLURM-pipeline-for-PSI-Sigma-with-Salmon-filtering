#!/bin/bash
#SBATCH --job-name=RNA_SEQ_PREP
#SBATCH --partition=glacier
#SBATCH --cpus-per-task=12
#SBATCH --mem=60G
#SBATCH --time=24:00:00
#SBATCH --output=/sci/labs/zvika.granot/segev.munitz/slurm_logs/%x_%j.out
#SBATCH --error=/sci/labs/zvika.granot/segev.munitz/slurm_logs/%x_%j.err

# -------------------------------------------------------------------
# HOW TO SUBMIT THIS JOB
#
# This is a single (non-array) preprocessing step.
#
# 1) From the directory that contains this .slurm file + the pipeline code:
#      sbatch sbatch_01_prep.slurm
#
# Expected inputs:
#   - samples.txt (one sample ID per line, no header)
#
# What it does:
#   - Runs: step1_prep_run.py
# -------------------------------------------------------------------

set -euo pipefail

# Ensure we run from the directory you submitted from
cd "$SLURM_SUBMIT_DIR"

python3 step1_prep_run.py
