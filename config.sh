# Shared paths. Load PLINK, GEMMA or ADMIXTURE in your environment before running.
set -euo pipefail
export COWPEA_ROOT="${COWPEA_ROOT:-$HOME/data/cowpealandrace_460}"
export COWPEA_OUTPUT="${COWPEA_OUTPUT:-$COWPEA_ROOT/reproduced}"
[[ "$COWPEA_ROOT" == /* && "$COWPEA_OUTPUT" == /* ]] || { echo 'Use absolute paths.' >&2; exit 1; }
data="$COWPEA_ROOT/manuscript_reanalysis"
geno="$data/genotype_rebuild_460_20261014/cowpea460_verified_geno005_maf001"
test -d "$data"
