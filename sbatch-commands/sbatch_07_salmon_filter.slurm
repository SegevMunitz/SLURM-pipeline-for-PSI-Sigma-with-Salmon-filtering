#!/bin/bash
#SBATCH --job-name=SALMON_FILTER
#SBATCH --partition=glacier
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=/sci/labs/zvika.granot/segev.munitz/slurm_logs/salmon_filter_%j.out
#SBATCH --error=/sci/labs/zvika.granot/segev.munitz/slurm_logs/salmon_filter_%j.err

# -------------------------------------------------------------------
# HOW TO SUBMIT THIS JOB
#
# This is a single (non-array) job.
# It should be run AFTER:
#   - Salmon quantification array (sbatch_06)
#   - PSI-Sigma run (sbatch_05)
#
# Submit:
#   sbatch sbatch_07_salmon_filter.slurm
#
# Optional chaining example (wait for Salmon array to finish):
#   N=$(grep -v '^\s*$' samples.txt | wc -l)
#   jid_salmon=$(sbatch --array=1-${N}%8 sbatch_06_salmon_array.slurm | awk '{print $4}')
#   sbatch --dependency=afterok:${jid_salmon} sbatch_07_salmon_filter.slurm
#
# What it does:
#   - Runs step7_psi_sigma_filtered_by_salmon.py
#   - Filters PSI-Sigma results using Salmon quantification (quant.sf)
# -------------------------------------------------------------------

set -euo pipefail

# Ensure we run from the submission directory
cd "$SLURM_SUBMIT_DIR"

python3 step7_psi_sigma_filtered_by_salmon.py
