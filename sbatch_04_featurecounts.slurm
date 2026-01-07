#!/bin/bash
#SBATCH --job-name=FEATURE_COUNTS
#SBATCH --partition=glacier
#SBATCH --cpus-per-task=12
#SBATCH --mem=60G
#SBATCH --time=24:00:00

############################ EDIT #############################

#SBATCH --output=/sci/labs/zvika.granot/segev.munitz/   /logs/%x_%j.out
#SBATCH --error=/sci/labs/zvika.granot/segev.munitz/    /logs/%x_%j.err

############################ END OF EDIT #############################

set -euo pipefail

OUTDIR="/abs/path/to/output_root"
GTF="/abs/path/to/annotation.gtf"

ls -1 "$OUTDIR"/bam/*/*.Aligned.sortedByCoord.out.bam >/dev/null

python step4_featurecounts.py \
    --out-root "$OUTDIR" \
    --annotation-gtf "$GTF" \
    --threads "$SLURM_CPUS_PER_TASK" \
    --stranded 0 \
    --paired-end
