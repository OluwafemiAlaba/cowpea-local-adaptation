#!/bin/bash -l
#SBATCH --job-name=cowpea460_verified_GWAS
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:20:00
#SBATCH --array=1-92%8
set -euo pipefail
module load gemma/0.98.5
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
project_root="${COWPEA_ROOT:-$HOME/data/cowpealandrace_460}"
root="$project_root/manuscript_reanalysis"
cd "$root/egwas_verified_460"
geno="$root/genotype_rebuild_460_20261014/cowpea460_verified_geno005_maf001"
sha256sum --status -c input/genotype_sha256.txt
sha256sum --status -c input/association_sha256.txt
column=${SLURM_ARRAY_TASK_ID:?Must run as a Slurm array task}
test "$column" -ge 1 && test "$column" -le 92
trait=$(sed -n "${column}p" input/trait_names.txt)
test -n "$trait"
test ! -e "results/${trait}.assoc.txt"
printf 'Trait: %s; phenotype column: %s\n' "$trait" "$column"
gemma -bfile "$geno" -p input/phenotypes_92.txt -n "$column" \
  -k kinship/cowpea460_centered.cXX.txt -maf 0.01 -miss 0.05 -hwe 0 \
  -r2 0.9999 -lmm 4 -silence -outdir results -o "$trait"
test -s "results/${trait}.assoc.txt"
grep -q 'number of analyzed individuals = 460' "results/${trait}.log.txt"
grep -q 'number of analyzed SNPs/var = 41654' "results/${trait}.log.txt"
test "$(wc -l < "results/${trait}.assoc.txt")" -eq 41655
sha256sum "results/${trait}.assoc.txt" "results/${trait}.log.txt" > "results/${trait}.sha256.txt"
