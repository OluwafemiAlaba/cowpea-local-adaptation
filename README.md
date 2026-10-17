# Cowpea local adaptation

Final cohort: 460 cowpea accessions; 41,654 array SNPs; 92 climate variables.

## Files

The GitHub archive contains 35 core analysis scripts and setup/validation helpers. The Dryad archive contains final processed inputs, the annotation needed for gene mapping, and key result tables. Extract both archives into the same directory; their `manuscript_reanalysis/` paths match.

Run `python3 verify_minimal.py` to check file hashes and core sample/marker counts. Run `python3 configure_paths.py` to preview historical HPCC path replacements, then use `--apply` in a working copy. Slurm partitions and module commands must match the execution environment. R package versions are recorded in the supplied session files. Python scripts require NumPy and pandas where imported; genome analyses require PLINK, GEMMA and ADMIXTURE; raster processing requires GDAL.

## Analysis order

1. Use `environment_460.csv` and `analysis_panel_460_metadata.csv` with the PLINK files in `genotype_rebuild_460_20261014`. Join by accession ID, never by row position. Southern monthly values are already shifted by six months: do not shift them again.
2. Structure: use `structure_LD02.prune.in`, `genetic_structure_verified_460/run_genetic_pca.R` and `admixture/` scripts. To rebuild the PCA input, export the pruned SNPs from the final PLINK files as VCF, keeping allele order; the PCA script takes VCF, two-column keep file, and output prefix. The supplied assignment table uses K=3 and Q >= 0.70.
3. Differentiation: run `genetic_differentiation_verified_460/calculate_fst.py` from the project root. Its `input/` directory expects links to the final PLINK files and `structure_LD02.prune.in` from the genotype folder; copies are omitted to avoid duplicating data.
4. RDA: run `run_rda.R`, `finish_rda.R` and `screen_candidates_3SD.R` from `rda_verified_460/`. Aligned inputs and the independent PLINK export are included. Outliers use +/-3 SD. Intermediate fitted RDS objects are regenerated.
5. eGWAS: in `egwas_verified_460/`, run kinship, the 92 association tasks, `009_post_gwas.sh`, and `parse_gwas.R`; then prepare tag inputs, calculate LD, select tags and map genes. Create `audit/`, `kinship/` and `results/` if absent. Full scan results are regenerated, not deposited. Per-trait and across-92 thresholds remain distinct; tags use r2 > 0.20 within 1 Mb on the same chromosome.
6. Compare candidates using `candidate_overlap_verified_460/` and build shared-gene tables using `shared_candidate_genes_460/`. Intermediate files needed by downstream table scripts are generated in the preceding steps.
7. GO: run `go_window20kb/map_candidates.py`; the approved foreground/background lists are under `topGO/combined/inputs/without_simple_*`. Run `run_without_simple.R` from that folder. It pools pRDA and per-trait eGWAS genes within +/-20 kb and excludes simple RDA. `prepare_union.py` documents gene-set construction but also generates the historical all-three comparison; use only the approved `without_simple` lists for the reported enrichment. topGO P values are nominal, not FDR-controlled discoveries.
8. Offset: in `GF_validation_design_review/`, run `run_comparison.R` with tasks 0-17, then the projection scripts. A small cached transformation and canonical raster-cell addresses are supplied for the site projection; full fitted forests are not bundled. Tasks 0, 6 and 12 fit full-data models; remaining tasks reproduce the completed geographic holdouts. Climate manifests specify source URLs, selected bands and the six-month future-data shift. Download/prepare climate rasters from those sources; they are intentionally not bundled. `prepare_landscape.R` also requires the baseline extraction produced by the original climate audit, which is supplied. Offset represents exploratory genomic mismatch, not demonstrated loss of fitness.

## Scope

This is a minimal analysis/data deposit, not an archive of every intermediate file or every figure layout. Saved models, full association scans, raster products, older datasets, manuscript copies and historical plot revisions remain in the private working archive. Core inputs and key results are supplied; a full clean-environment model rerun has not been performed during packaging. Scripts retain recorded analysis assumptions and may require environment/path setup. The original broad package is superseded for deposition by this minimal package.

Original genotype resource: https://doi.org/10.6086/D19Q37. Passport/environment records were supplied through IITA. Retain attribution and confirm redistribution terms for cowpea annotation, GO mappings and climate data before public release. Private code repository: https://github.com/OluwafemiAlaba/cowpea-local-adaptation. The Dryad data archive remains local; code license and Dryad DOI are pending.

Corresponding author: Daniel Koenig, dkoenig@ucr.edu.
The authors declare no competing interests.
