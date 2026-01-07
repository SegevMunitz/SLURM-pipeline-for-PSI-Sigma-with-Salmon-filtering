#!/bin/bash
#SBATCH --job-name=STAR_ALIGN
#SBATCH --partition=glacier
#SBATCH --cpus-per-task=12
#SBATCH --mem=60G
#SBATCH --time=24:00:00
#SBATCH --array=1-1%6
############################ EDIT #############################

#SBATCH --output=/sci/labs/zvika.granot/segev.munitz/   /logs/%x_%j.out
#SBATCH --error=/sci/labs/zvika.granot/segev.munitz/    /logs/%x_%j.err

############################ END OF EDIT #############################

set -euo pipefail

OUTDIR="/abs/path/to/output_root"             # must match paths.OUTDIR
INDEX="/abs/path/to/star_index/mm10_rl100"    # paths.STAR_INDEX_DIR

mkdir -p "$OUTDIR/logs"

python step2_star_align_core.py \
  --genome-dir "$INDEX" \
  --out-dir "$OUTDIR" \
  --samples "$OUTDIR/samples.txt" \
  --threads "$SLURM_CPUS_PER_TASK"


# submit with:

# OUTDIR="/abs/path/to/output_root"          # must match paths.OUTDIR
# N=$(wc -l < "$OUTDIR/samples.txt")

# jid=$(sbatch --array=1-${N}%6 sbatch_02_star_array.slurm | awk '{print $4}')
# sbatch --dependency=afterok:${jid} sbatch_02b_make_groups.slurm
