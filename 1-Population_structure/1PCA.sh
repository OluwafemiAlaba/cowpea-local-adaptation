#!/bin/bash
# Export 4,194 LD-pruned SNPs for genetic PCA.
source ./config.sh
out="$COWPEA_OUTPUT/structure"
mkdir -p "$out"
plink --bfile "$geno" --extract "$data/genotype_rebuild_460_20261014/structure_LD02.prune.in" \
  --allow-extra-chr --keep-allele-order --recode vcf-iid --out "$out/primary460"
awk '{print $1, $2}' "$geno.fam" > "$out/keep.txt"
Rscript 1-Population_structure/2PCA.R "$out/primary460.vcf" "$out/keep.txt" "$out/primary460"
