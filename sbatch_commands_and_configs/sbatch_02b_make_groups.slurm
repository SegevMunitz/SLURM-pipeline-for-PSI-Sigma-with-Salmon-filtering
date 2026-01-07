#SBATCH --job-name=MAKE_GROUPS
#SBATCH --partition=glacier
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:10:00
############################ EDIT #############################

#SBATCH --output=/sci/labs/zvika.granot/segev.munitz/   /logs/%x_%j.out
#SBATCH --error=/sci/labs/zvika.granot/segev.munitz/    /logs/%x_%j.err

############################ END OF EDIT #############################

set -euo pipefail

OUTDIR="/abs/path/to/output_root"

python step2b_make_group_files.py \
  --out-dir "$OUTDIR" \
  --samples "$OUTDIR/samples.txt" \
  --map "$OUTDIR/groups.tsv"
