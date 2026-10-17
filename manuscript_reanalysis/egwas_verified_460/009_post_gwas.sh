#!/bin/bash
#SBATCH --job-name=cowpea460_postgwas
#SBATCH --partition=short
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:10:00
# Adapted from the archived 009_post_gwas.sh; original author unspecified.
# Local: bash 009_post_gwas.sh [absolute egwas_verified_460 directory]
# HPCC: sbatch 009_post_gwas.sh /absolute/path/to/egwas_verified_460
# This only summarizes existing scans. No association or LD model is fitted.
set -euo pipefail
export LC_ALL=C
root=${1:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
results="$root/results"
out="$root/post_gwas"
traits="$root/input/trait_names.txt"
test -s "$traits"
mkdir -p "$out"
stage=$(mktemp -d "$out/.stage.XXXXXX")
trap 'rm -rf "$stage"' EXIT
lowp="$stage/cowpea460_gwas_snps_lowp.txt"
tested="$stage/cowpea460_gwas_snps_tested.txt"
pve="$stage/pve.txt"
: > "$lowp"
: > "$tested"
: > "$pve"
trait_count=0
while IFS= read -r trait || [[ -n "$trait" ]]; do
    [[ "$trait" =~ ^[A-Za-z0-9_]+$ ]] || { printf 'Invalid trait name\n' >&2; exit 1; }
    test -s "$results/$trait.assoc.txt"
    test -s "$results/$trait.log.txt"
    # Preserve all original GEMMA columns and prefix the trait, as in the old
    # lowp file. Use the named LRT column and skip the header explicitly.
    awk -v trait="$trait" -v lowp="$lowp" -v tested="$tested" '
        BEGIN { FS=OFS="\t"; expected="chr rs ps n_miss allele1 allele0 af beta se logl_H1 l_remle l_mle p_wald p_lrt p_score"; split(expected,names," ") }
        NR==1 {
            if (NF!=15) {bad=1; exit 1}
            for(i=1;i<=15;i++) {
                if($i!=names[i]) {bad=1; exit 1}
                if($i=="p_lrt") pcol=i
            }
            next
        }
        NF {
            if (NF!=15 || $pcol !~ /^[0-9.eE+-]+$/ || $pcol+0<=0 || $pcol+0>1 || seen[$2]++) {bad=1; exit 1}
            n++
            if ($pcol+0 < 0.0001) print trait,$0 >> lowp
        }
        END {
            if(bad || n!=41654 || !pcol) exit 1
            print n,trait >> tested
        }
    ' "$results/$trait.assoc.txt" || { printf 'Invalid association records: %s\n' "$trait" >&2; exit 1; }
    awk -v trait="$trait" '
        /## pve estimate in the null model =/ {
            value=$NF; n++
            if(value !~ /^[0-9.eE+-]+$/ || value+0<0 || value+0>1) bad=1
        }
        END {if(n!=1 || bad) exit 1; printf "%s\t%s\n",trait,value}
    ' "$results/$trait.log.txt" >> "$pve"
    trait_count=$((trait_count+1))
done < "$traits"
test "$trait_count" -eq 92
test "$(cut -f2 "$tested" | sort -u | wc -l | tr -d ' ')" -eq 92
for name in cowpea460_gwas_snps_lowp.txt cowpea460_gwas_snps_tested.txt pve.txt; do
    mv "$stage/$name" "$out/$name"
done
printf 'Completed %s traits; 41,654 SNPs per trait; %s low-P SNP-trait rows (LRT P < 0.0001).\n' \
    "$trait_count" "$(wc -l < "$out/cowpea460_gwas_snps_lowp.txt" | tr -d ' ')"
