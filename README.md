# Cowpea local adaptation and genomic offset

R and shell scripts for studying local adaptation and genomic offset in 460 cowpea accessions using 41,654 array SNPs and 92 climate variables.

## Requirements

R packages: adegenet, ade4, pegas, vegan, adespatial, spdep, geosphere, data.table, jsonlite, topGO and gradientForest.

Software: PLINK 1.9, GEMMA 0.98.5 and ADMIXTURE 1.3.0.

## Data and paths

Extract the Dryad data into `~/data/cowpealandrace_460`, keeping the `manuscript_reanalysis/` folder. Run the commands below from this repository. Results are saved in `~/data/cowpealandrace_460/reproduced/`.

To use other locations, set `COWPEA_ROOT` and `COWPEA_OUTPUT` to absolute paths. Set `COWPEA_USE_DEPOSITED=1` to use available saved results for later steps.

## Using your own data

Update the paths in `config.R` and `config.sh`, then change input filenames, sample and SNP counts, predictor lists and analysis settings in the relevant scripts. Keep accession and SNP identifiers consistent across input files. Counts and settings in these scripts refer to the cowpea dataset.

## 1. Population structure

PCA uses 4,194 SNPs pruned for linkage disequilibrium. ADMIXTURE runs K = 1–10 with ten replicates and fivefold cross-validation. Cluster assignments use Q ≥ 0.70 at K = 3.

```sh
bash 1-Population_structure/1PCA.sh
bash 1-Population_structure/3ADMIXTURE.sh prepare
for task in $(seq 0 99); do bash 1-Population_structure/3ADMIXTURE.sh "$task"; done
```

## 2. Redundancy analysis

Simple RDA, partial RDA and variance partitioning use 20,027 pruned SNPs. Partial RDA accounts for genetic PC1–PC4 and spatial eigenvectors. Candidate SNPs are identified at ±3 SD on supported axes.

```sh
Rscript 2-RDA/1RDA.R
Rscript 2-RDA/2Outliers.R
```

## 3. Environmental GWAS

GEMMA tests 92 climate variables using mixed models. Likelihood-ratio P values are compared with Bonferroni thresholds of 0.05/41,654 per trait and 0.05/(41,654 × 92) across all traits. Tag SNPs represent directly linked significant SNPs on the same chromosome, less than 1 Mb apart, with r² > 0.20.

```sh
bash 3-eGWAS/1GEMMA.sh kinship
for trait in $(seq 1 92); do bash 3-eGWAS/1GEMMA.sh "$trait"; done
Rscript 3-eGWAS/2Parse.R
bash 3-eGWAS/3LD.sh
Rscript 3-eGWAS/4Tags.R
```

## 4. Candidate genes and GO enrichment

Genes within ±20 kb of pRDA and eGWAS candidate SNPs are combined. The background uses the same window around all tested SNPs. Biological Process enrichment is tested with topGO using the elim Fisher test.

```sh
Rscript 4-Gene_enrichment/1Candidate_genes.R
Rscript 4-Gene_enrichment/2topGO.R
```

## 5. Genomic offset

Gradient Forest uses six climate predictors and three SNP sets: pRDA, pRDA plus eGWAS, and simple RDA plus pRDA plus eGWAS. Genomic offset is calculated at 348 collection locations across three climate models, three SSPs and three periods.

```sh
Rscript 5-Genomic_offset/1GradientForest.R
Rscript 5-Genomic_offset/2Genomic_offset.R
```

Monthly climate values for southern accessions have already been shifted by six months. Do not repeat this shift. Bioclimatic variables are unchanged.

## References

- [Cowpea genotype data](https://doi.org/10.6086/D19Q37)
- [Fiscus et al. (2024)](https://doi.org/10.1093/g3journal/jkae071)

## Contact

Oluwafemi A. Alaba: [oluwafea@ucr.edu](mailto:oluwafea@ucr.edu)
