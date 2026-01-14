#!/bin/bash
#SBATCH --job-name=SALMON
#SBATCH --partition=glacier
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --array=1-1%8
#SBATCH --output=/sci/labs/zvika.granot/segev.munitz/slurm_logs/salmon_%A_%a.out
#SBATCH --error=/sci/labs/zvika.granot/segev.munitz/slurm_logs/salmon_%A_%a.err

# -------------------------------------------------------------------
# HOW TO SUBMIT THIS JOB (ARRAY)
#
# This is an array job: one task per sample in samples.txt
#
# 1) Make sure samples file exists (one sample ID per line, no header):
#      samples.txt
#
# 2) Count samples:
#      N=$(grep -v '^\s*$' samples.txt | wc -l)
#
# 3) Submit (recommended throttled):
#      sbatch --array=1-${N}%8 sbatch_06_salmon_array.slurm
#
# Notes:
# - SLURM_ARRAY_TASK_ID runs from 1..N
#
# Important:
# - Task 1 will build the Salmon index if it does not already exist,
#   then proceed to quantification.
# -------------------------------------------------------------------

set -euo pipefail

# Ensure we run from the submission directory
cd "$SLURM_SUBMIT_DIR"

python3 step6_salmon_quant_core.py
