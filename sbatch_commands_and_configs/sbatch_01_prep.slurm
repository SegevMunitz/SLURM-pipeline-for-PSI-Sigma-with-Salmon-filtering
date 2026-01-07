#!/bin/bash
#SBATCH --job-name=RNA_SEQ_PREP
#SBATCH --partition=glacier
#SBATCH --cpus-per-task=12
#SBATCH --mem=60G
#SBATCH --time=24:00:00
############################ EDIT #############################

#SBATCH --output=/sci/labs/zvika.granot/segev.munitz/   /logs/%x_%j.out
#SBATCH --error=/sci/labs/zvika.granot/segev.munitz/    /logs/%x_%j.err

############################ END OF EDIT #############################

set -euo pipefail

python step1_prep_run.py
