#!/usr/bin/env bash

################# HOW TO SUBMIT ALL JOBS IN ORDER #################

# This script will submit all the sbatch jobs in the correct order,
# with appropriate dependencies.
# It assumes you have a samples.txt file in the current directory.
#
# To use:
# 1) Save this script as sbatch_commands/submit_all.sh
# 2) Make it executable and run:

#           chmod +x submit_all.sh
#           ./submit_all.sh

###################################################################

set -euo pipefail

# ---- config ----
MAXP_STAR=8
MAXP_SALMON=8

# Ensure log directory exists
LOGDIR="/sci/labs/zvika.granot/segev.munitz/slurm_logs"
mkdir -p "$LOGDIR" 

# Count samples
N=$(grep -v '^\s*$' samples.txt | wc -l)
if [[ "$N" -lt 1 ]]; then
  echo "ERROR: samples.txt has no samples"
  exit 1
fi
echo "Samples: N=$N"

# 1) Prep
jid_prep=$(sbatch sbatch_01_prep.slurm | awk '{print $4}')
echo "prep: $jid_prep"

# 2) STAR array after prep
jid_star=$(sbatch --dependency=afterok:${jid_prep} --array=1-${N}%${MAXP_STAR} sbatch_02_star_array.slurm | awk '{print $4}')
echo "star: $jid_star"

# 3) Make groups after STAR finishes
jid_groups=$(sbatch --dependency=afterok:${jid_star} sbatch_03_make_groups.slurm | awk '{print $4}')
echo "groups: $jid_groups"

# 4) featureCounts after STAR finishes (or after groups; your choice)
jid_fc=$(sbatch --dependency=afterok:${jid_star} sbatch_04_featurecounts.slurm | awk '{print $4}')
echo "featureCounts: $jid_fc"

# 5) PSI-Sigma after groups + featureCounts
jid_psisigma=$(sbatch --dependency=afterok:${jid_groups}:${jid_fc} sbatch_05_run_psi_sigma.slurm | awk '{print $4}')
echo "psi-sigma: $jid_psisigma"

# 6) Salmon array after STAR (needs BAMs)
jid_salmon=$(sbatch --dependency=afterok:${jid_star} --array=1-${N}%${MAXP_SALMON} sbatch_06_salmon_array.slurm | awk '{print $4}')
echo "salmon: $jid_salmon"

# 7) Salmon-filter after PSI-Sigma + Salmon
jid_filter=$(sbatch --dependency=afterok:${jid_psisigma}:${jid_salmon} sbatch_07_salmon_filter.slurm | awk '{print $4}')
echo "salmon_filter: $jid_filter"

echo ""
echo "Submitted all jobs."
echo "Final job (filter) will run after everything succeeds: $jid_filter"
