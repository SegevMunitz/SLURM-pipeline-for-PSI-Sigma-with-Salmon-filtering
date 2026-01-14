#!/bin/bash
#SBATCH --job-name=STAR_ALIGN
#SBATCH --partition=glacier
#SBATCH --cpus-per-task=12
#SBATCH --mem=60G
#SBATCH --time=24:00:00
#SBATCH --array=1-1%8
#SBATCH --output=/sci/labs/zvika.granot/segev.munitz/slurm_logs/%x_%A_%a.out
#SBATCH --error=/sci/labs/zvika.granot/segev.munitz/slurm_logs/%x_%A_%a.err

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
#      sbatch --array=1-${N}%8 sbatch_02_star_array.slurm
#
# Notes:
# - SLURM_ARRAY_TASK_ID runs from 1..N
#
# Optional chaining example (run group-making after STAR finishes):
#   jid=$(sbatch --array=1-${N}%8 sbatch_02_star_array.slurm | awk '{print $4}')
#   sbatch --dependency=afterok:${jid} sbatch_02b_make_groups.slurm
# -------------------------------------------------------------------

set -euo pipefail

# Ensure we run from the directory you submitted from
cd "$SLURM_SUBMIT_DIR"

python3 step2_star_align_core.py
