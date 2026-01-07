#!/bin/bash
#SBATCH --job-name=PSI_SIGMA_RUN
#SBATCH --partition=glacier
#SBATCH --cpus-per-task=16
#SBATCH --mem=60G
#SBATCH --time=24:00:00

############################ EDIT #############################

#SBATCH --output=/sci/labs/zvika.granot/segev.munitz/   /logs/%x_%j.out
#SBATCH --error=/sci/labs/zvika.granot/segev.munitz/    /logs/%x_%j.err

############################ END OF EDIT #############################

set -euo pipefail

echo "=== $(date) | job=$SLURM_JOB_ID node=$(hostname) ==="
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo ""

# Optional: load modules here (cluster-specific)
# module load STAR samtools perl python

# Run PSI-Sigma wrapper
python step5_run_psi_sigma.py
