#!/bin/bash -l
#SBATCH --job-name=cowpea460_verified_K
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:15:00
set -euo pipefail
module load gemma/0.98.5
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
project_root="${COWPEA_ROOT:-$HOME/data/cowpealandrace_460}"
root="$project_root/manuscript_reanalysis"
cd "$root/egwas_verified_460"
geno="$root/genotype_rebuild_460_20261014/cowpea460_verified_geno005_maf001"
sha256sum -c input/genotype_sha256.txt
sha256sum -c input/kinship_input_sha256.txt
test ! -e kinship/cowpea460_centered.cXX.txt
mkdir -p kinship logs
gemma -h 2 > logs/gemma_file_options.txt
gemma -h 3 > logs/gemma_QC_options.txt
gemma -bfile "$geno" -p input/climate_92_raw_units.txt -n 1 \
  -maf 0.01 -miss 0.05 -gk 1 -outdir kinship -o cowpea460_centered
test -s kinship/cowpea460_centered.cXX.txt
grep -q 'number of analyzed individuals = 460' kinship/cowpea460_centered.log.txt
grep -q 'number of analyzed SNPs/var = 41654' kinship/cowpea460_centered.log.txt
sha256sum kinship/cowpea460_centered.cXX.txt > kinship/output_sha256.txt
