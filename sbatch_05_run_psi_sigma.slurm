#!/bin/bash
#SBATCH --job-name=PSI_SIGMA_RUN
#SBATCH --partition=glacier
#SBATCH --cpus-per-task=16
#SBATCH --mem=60G
#SBATCH --time=24:00:00
#SBATCH --output=/sci/labs/zvika.granot/segev.munitz/slurm_logs/%x_%j.out
#SBATCH --error=/sci/labs/zvika.granot/segev.munitz/slurm_logs/%x_%j.err

# -------------------------------------------------------------------
# HOW TO SUBMIT THIS JOB
#
# This is a single (non-array) job.
# It should be run AFTER:
#   - STAR alignment (sbatch_02)
#   - Group file generation (sbatch_03)
#   - featureCounts (sbatch_04)
#
# Submit:
#   sbatch sbatch_05_run_psi_sigma.slurm
#
# Optional chaining example:
#   jid=$(sbatch sbatch_04_featurecounts.slurm | awk '{print $4}')
#   sbatch --dependency=afterok:${jid} sbatch_05_run_psi_sigma.slurm
#
# What it does:
#   - Runs step5_run_psi_sigma.py
#   - Executes PSI-Sigma on prepared group files
# -------------------------------------------------------------------

set -euo pipefail

# Ensure we run from the submission directory
cd "$SLURM_SUBMIT_DIR"

echo "=== $(date) | job=$SLURM_JOB_ID node=$(hostname) ==="
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo ""

# Run PSI-Sigma wrapper
python3 step5_run_psi_sigma.py
