#!/bin/bash
# LD against all tested markers; additive export is only for exact tag-pair checks.
source ./config.sh
base="$COWPEA_OUTPUT/egwas/post_gwas"
snps="$base/parsed/per_trait/cowpea460_gwas_hitlst_uniq_snps.txt"
if [[ ! -f "$snps" && "${COWPEA_USE_DEPOSITED:-0}" == 1 ]]; then
  snps="$data/egwas_verified_460/post_gwas/parsed/per_trait/cowpea460_gwas_hitlst_uniq_snps.txt"
fi
test -s "$snps"
mkdir -p "$base/tags"
plink --bfile "$geno" --allow-extra-chr --threads 1 --memory 1536 \
  --ld-snp-list "$snps" --r2 gz --ld-window 41654 --ld-window-kb 1000 --ld-window-r2 0 \
  --out "$base/tags/cowpea460_candidate_local_LD"
plink --bfile "$geno" --allow-extra-chr --keep-allele-order --extract "$snps" \
  --recode A --out "$base/tags/candidate_dosages"
