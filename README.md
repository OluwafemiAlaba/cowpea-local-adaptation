# Cowpea local adaptation and genomic offset

Selected R and shell scripts supporting the cowpea manuscript: **460 accessions, 41,654 array SNPs and 92 climate variables**. The numbered layout follows the organization of the [Populus repository](https://github.com/jingwanglab/Populus_genomic_prediction_climate_vulnerability); the analyses and settings are from this cowpea study.

## Setup

Extract the associated Dryad data into `~/data/cowpealandrace_460`, retaining `manuscript_reanalysis/`. Run commands from this repository's root. Set `COWPEA_ROOT` and `COWPEA_OUTPUT` to absolute paths to change the data and output locations. By default, new results go into `reproduced/` beside the deposited inputs.

Requires R with adegenet, ade4, pegas, vegan, adespatial, spdep, geosphere, data.table, jsonlite, topGO and gradientForest; PLINK 1.9; GEMMA 0.98.5; and ADMIXTURE 1.3.0. Recorded R package versions accompany the data.

## 1. Population structure

- `1-Population_structure/1PCA.sh` exports the 4,194 LD-pruned markers and calls `2PCA.R`: mean replacement of missing values, allele-frequency scaling and ade4 PCA.
- `3ADMIXTURE.sh`: K = 1–10, ten replicates and fivefold cross-validation. Final K = 3 assignments at Q ≥ 0.70 are deposited; run selection and label alignment scripts are omitted.

```sh
bash 1-Population_structure/1PCA.sh
bash 1-Population_structure/3ADMIXTURE.sh prepare
for task in $(seq 0 99); do bash 1-Population_structure/3ADMIXTURE.sh "$task"; done
```

## 2. Redundancy analysis

`1RDA.R` uses 20,027 pruned SNPs, mode imputation and SNP standardization. It screens climate predictors, selects spatial MEMs, fits simple RDA and pRDA, and partitions variation. pRDA conditions on genetic PC1–PC4 and selected MEMs. `2Outliers.R` screens supported axes at **±3 SD** and assigns putative candidates to their strongest climate correlation.

```sh
Rscript 2-RDA/1RDA.R
Rscript 2-RDA/2Outliers.R
```

## 3. Environmental GWAS

`1GEMMA.sh` calculates centered kinship and runs 92 mixed models. `2Parse.R` uses likelihood-ratio P values at **0.05/41,654 per trait** and **0.05/(41,654 × 92) across all scans**. The P < 0.0001 extract is only an intermediate.

`3LD.sh` calculates LD. `4Tags.R` groups each lowest-P SNP with its directly linked significant neighbours: same chromosome, **distance <1 Mb and r² >0.20**. Singletons are retained; LD is checked against observed genotypes.

```sh
bash 3-eGWAS/1GEMMA.sh kinship
for trait in $(seq 1 92); do bash 3-eGWAS/1GEMMA.sh "$trait"; done
Rscript 3-eGWAS/2Parse.R
bash 3-eGWAS/3LD.sh
Rscript 3-eGWAS/4Tags.R
```

The ADMIXTURE and GEMMA loops can instead run as Slurm arrays, with indices 0–99 and 1–92, respectively.

## 4. Candidate genes and GO enrichment

`1Candidate_genes.R` pools pRDA and per-trait eGWAS genes within **±20 kb** of candidate SNPs, excluding simple RDA. Background genes use the same window around all tested array SNPs. Genes are deduplicated. `2topGO.R` tests Biological Process terms using elim Fisher tests; nominal P values are exploratory, and adjusted values are retained separately.

```sh
Rscript 4-Gene_enrichment/1Candidate_genes.R
Rscript 4-Gene_enrichment/2topGO.R
```

## 5. Genomic offset

`1GradientForest.R` fits the three recorded candidate sets using 500 trees, mtry = 2 and seeds 460, 1460 and 2460. The six predictors are prec_08, prec_01, bio_15, bio_3, tmin_10 and bio_19. `2Genomic_offset.R` calculates Euclidean distances between current and future transformed climates for 348 collection locations, three GCMs, three SSPs and three periods, and flags predictors outside sampled ranges.

```sh
Rscript 5-Genomic_offset/1GradientForest.R
Rscript 5-Genomic_offset/2Genomic_offset.R
```

The supplied future tables are already extracted and aligned. Southern monthly values use a six-month shift (January corresponds to July); **do not shift these tables again**. Bioclimatic variables are not shifted. Offset is an exploratory measure of genomic mismatch, not a measured loss of fitness.

## Scope and checks

Set `COWPEA_USE_DEPOSITED=1` to allow available deposited intermediates when new results are absent. RDA outlier screening requires the state saved by `1RDA.R`; topGO requires lists from `1Candidate_genes.R`. Raster downloads, figure rendering and optional checks are omitted.

The simplified candidate screening, eGWAS parsing, tag grouping, GO gene mapping and site-offset calculations were compared with saved results. A complete refit of all models was not performed for this code revision. Original scripts remain in Git history.

Genotype resource: [UCR cowpea data](https://doi.org/10.6086/D19Q37). Tag grouping follows the stated cowpea LD criterion informed by [Fiscus et al. (2024)](https://doi.org/10.1093/g3journal/jkae071). GO code retains the topGO approach adapted from Christopher Fiscus's original script.
