#!/bin/bash -l
#SBATCH --job-name=cowpea460_tagLD
#SBATCH --partition=short
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:15:00
set -euo pipefail
module load plink
root=${1:?Supply absolute egwas_verified_460 directory}
geno="$root/../genotype_rebuild_460_20261014/cowpea460_verified_geno005_maf001"
out="$root/post_gwas/tags"
test "$(wc -l < "$geno.fam")" -eq 460
awk '$3!=0 || $4!=0 {exit 1}' "$geno.fam"
test "$(wc -l < "$out/per_trait_unique_snps.txt")" -gt 0
sha256sum "$geno.bed" "$geno.bim" "$geno.fam" > "$out/genotype_source.sha256"
# Do not --extract the hit list: measure each significant focal SNP against
# all typed QC markers in its local window. Retain all r2 for transparent
# downstream selection, rather than relying on the default reporting cutoff.
plink --bfile "$geno" --allow-extra-chr --threads 1 --memory 1536 \
  --ld-snp-list "$out/per_trait_unique_snps.txt" --r2 gz \
  --ld-window 41654 --ld-window-kb 1000 --ld-window-r2 0 \
  --out "$out/cowpea460_candidate_local_LD"
