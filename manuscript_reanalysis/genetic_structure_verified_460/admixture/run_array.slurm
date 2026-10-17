#!/bin/bash -l
#SBATCH --job-name=cowpea_admix460
#SBATCH --partition=short
#SBATCH --cpus-per-task=2
#SBATCH --mem=2G
#SBATCH --time=01:00:00
#SBATCH --array=0-99%10
set -euo pipefail
module load admixture/1.3.0
project_root="${COWPEA_ROOT:-$HOME/data/cowpealandrace_460}"
base="$project_root/manuscript_reanalysis/genetic_structure_verified_460/admixture"
k=$((SLURM_ARRAY_TASK_ID / 10 + 1))
rep=$((SLURM_ARRAY_TASK_ID % 10 + 1))
seed=$((460000 + k * 100 + rep))
work="$base/K${k}/rep${rep}"
mkdir -p "$work"
cd "$work"
for ext in bed bim fam; do ln -s "$base/input/primary460.$ext" "primary460.$ext"; done
printf 'K=%s\nreplicate=%s\nseed=%s\nCV_folds=5\n' "$k" "$rep" "$seed" > run_parameters.txt
admixture --cv=5 -j2 -s "$seed" primary460.bed "$k" > admixture.log 2>&1
[[ -s "primary460.${k}.Q" && -s "primary460.${k}.P" ]]
printf 'COMPLETE\n' > COMPLETE
