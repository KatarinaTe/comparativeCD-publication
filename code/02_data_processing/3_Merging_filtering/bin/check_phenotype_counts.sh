#!/usr/bin/env bash
# check_phenotype_counts.sh — cross-phenotype dog/SNP-count check
# Source: 02_data_processing/3_Merging_filtering/02_create_ALLFAM_check_stats.sh, lines 26-84
# ("check SNPs and dogs in each GWAS-dataset" terminal section)
#
# The original ran `wc -l` by hand on each of the 17 (3 factor + 14 item) QC6 .fam/.bim files and
# transcribed the results into a comment-block table. Reproduced here as the same `wc -l` calls
# followed by a computed table, rather than the original's hardcoded copy-pasted numbers.
#
# Args: <fam1> [fam2 ...] -- <bim1> [bim2 ...] -- <id1> [id2 ...]
# (fam/bim/id lists must be the same length and in the same phenotype order)

set -euo pipefail

fams=()
bims=()
ids=()
section=0
for arg in "$@"; do
    if [ "$arg" = "--" ]; then
        section=$((section + 1))
        continue
    fi
    case $section in
        0) fams+=("$arg") ;;
        1) bims+=("$arg") ;;
        2) ids+=("$arg") ;;
    esac
done

echo "# N dogs:"
for f in "${fams[@]}"; do
    wc -l "$f"
done

echo
echo "# N SNPs:"
for b in "${bims[@]}"; do
    wc -l "$b"
done

echo
printf "%s\t%s\t%s\n" "QC6" "Dogs" "SNPs"
for i in "${!ids[@]}"; do
    ndogs=$(wc -l < "${fams[$i]}")
    nsnps=$(wc -l < "${bims[$i]}")
    printf "%s\t%s\t%s\n" "${ids[$i]}" "$ndogs" "$nsnps"
done
