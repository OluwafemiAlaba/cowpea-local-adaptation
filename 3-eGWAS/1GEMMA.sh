#!/bin/bash
# bash 3-eGWAS/1GEMMA.sh kinship; then bash 3-eGWAS/1GEMMA.sh TRAIT_INDEX (1..92).
# Or: sbatch --array=1-92 3-eGWAS/1GEMMA.sh (choose resources for your cluster).
source ./config.sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
stage="$data/egwas_verified_460"
out="$COWPEA_OUTPUT/egwas"
mkdir -p "$out/kinship" "$out/results"
mode="${1:-${SLURM_ARRAY_TASK_ID:-}}"
if [[ "$mode" == kinship ]]; then
  test ! -e "$out/kinship/cowpea460_centered.cXX.txt"
  # The first complete phenotype column supplies the same 460-sample inclusion mask.
  gemma -bfile "$geno" -p "$stage/input/phenotypes_92.txt" -n 1 \
    -maf 0.01 -miss 0.05 -gk 1 -outdir "$out/kinship" -o cowpea460_centered
  log="$out/kinship/cowpea460_centered.log.txt"
else
  [[ "$mode" =~ ^[0-9]+$ ]] && (( mode >= 1 && mode <= 92 )) || { echo 'Use kinship or a trait index 1..92.' >&2; exit 1; }
  trait=$(sed -n "${mode}p" "$stage/input/trait_names.txt")
  test -n "$trait" && test ! -e "$out/results/${trait}.assoc.txt"
  gemma -bfile "$geno" -p "$stage/input/phenotypes_92.txt" -n "$mode" \
    -k "$out/kinship/cowpea460_centered.cXX.txt" -maf 0.01 -miss 0.05 -hwe 0 \
    -r2 0.9999 -lmm 4 -silence -outdir "$out/results" -o "$trait"
  test "$(wc -l < "$out/results/${trait}.assoc.txt")" -eq 41655
  log="$out/results/${trait}.log.txt"
fi
grep -q 'number of analyzed individuals = 460' "$log"
grep -q 'number of analyzed SNPs/var = 41654' "$log"
