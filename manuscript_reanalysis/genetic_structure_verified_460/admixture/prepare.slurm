#!/bin/bash -l
#SBATCH --job-name=cowpea_admix_prep
#SBATCH --partition=short
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:10:00
set -euo pipefail
module load plink/1.90b6.25
project_root="${COWPEA_ROOT:-$HOME/data/cowpealandrace_460}"
base="$project_root/manuscript_reanalysis"
cd "$base/genetic_structure_verified_460/admixture"
mkdir -p input
plink --bfile "$base/genotype_rebuild_460_20261014/cowpea460_verified_geno005_maf001" --extract "$base/genotype_rebuild_460_20261014/structure_LD02.prune.in" --allow-extra-chr --keep-allele-order --make-bed --memory 1500 --threads 1 --out input/primary460
# ADMIXTURE requires integer chromosome names. Only the BIM chromosome labels change.
cp input/primary460.bim input/primary460.Vu_labels.bim
awk 'BEGIN{OFS="\t"} {if($1 !~ /^Vu(0[1-9]|1[01])$/)exit 1; sub(/^Vu/,"",$1);$1=$1+0;print}' input/primary460.Vu_labels.bim > input/primary460.bim
[[ $(wc -l < input/primary460.fam) -eq 460 ]]
[[ $(wc -l < input/primary460.bim) -eq 4194 ]]
sha256sum input/primary460.* > input_sha256.txt
