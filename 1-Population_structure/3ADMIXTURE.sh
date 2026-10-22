#!/bin/bash
# First: bash 1-Population_structure/3ADMIXTURE.sh prepare
# Then: bash 1-Population_structure/3ADMIXTURE.sh TASK, for TASK 0..99.
# Alternatively submit with sbatch --array=0-99 (resources are site-specific).
source ./config.sh
base="$COWPEA_OUTPUT/structure/admixture"
task="${1:-${SLURM_ARRAY_TASK_ID:-}}"
if [[ "$task" == prepare ]]; then
  mkdir -p "$base/input"
  plink --bfile "$geno" --extract "$data/genotype_rebuild_460_20261014/structure_LD02.prune.in" \
    --allow-extra-chr --keep-allele-order --make-bed --out "$base/input/primary460"
  cp "$base/input/primary460.bim" "$base/input/primary460.Vu_labels.bim"
  awk 'BEGIN{OFS="\t"} {if($1 !~ /^Vu(0[1-9]|1[01])$/)exit 1;sub(/^Vu/,"",$1);$1=$1+0;print}' \
    "$base/input/primary460.Vu_labels.bim" > "$base/input/primary460.bim"
  test "$(wc -l < "$base/input/primary460.bim")" -eq 4194
  test "$(wc -l < "$base/input/primary460.fam")" -eq 460
  exit 0
fi
[[ "$task" =~ ^[0-9]+$ ]] && (( task >= 0 && task < 100 )) || { echo 'Use prepare or a task index 0..99.' >&2; exit 1; }
k=$((task / 10 + 1)); rep=$((task % 10 + 1)); seed=$((460000 + k * 100 + rep))
work="$base/K${k}/rep${rep}"
test ! -e "$work/COMPLETE"
mkdir -p "$work"
cd "$work"
for ext in bed bim fam; do ln -sf "$base/input/primary460.$ext" "primary460.$ext"; done
printf 'K=%s\nreplicate=%s\nseed=%s\nCV_folds=5\n' "$k" "$rep" "$seed" > run_parameters.txt
admixture --cv=5 -j2 -s "$seed" primary460.bed "$k" > admixture.log 2>&1
test -s "primary460.${k}.Q" && test -s "primary460.${k}.P"
printf 'COMPLETE\n' > COMPLETE
