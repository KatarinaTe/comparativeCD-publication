#!/usr/bin/env bash
# reformat_polmm_output.sh — splits GRAB.Marker's colon-packed marker ID columns into separate
# CHR/Position/Chr/bp/A1/A2 fields, passing every other column through unchanged.
# Source: 04_gwas/2_polmm/POLMMgrab_PER_ITEM.sh, line 44 ("polmm sumstats was modified like this")
#
# GRAB.Marker's own output column 1 is "chr:pos"-formatted, column 2 is "chr:pos:A1:A2"-formatted
# — this reproduces the original's exact awk one-liner splitting both, unchanged otherwise.
#
# Column 1 is just an echo of the input .bim file's marker ID (GRAB doesn't compute it) — the
# "chr:pos" assumption only holds because both upstream genotype pipelines
# (01_mapping/2_Axiom_imputation) set variant IDs via `bcftools annotate --set-id '%CHROM:%POS'`
# before merging. If that upstream convention were ever missing, this awk step would silently
# produce an empty Position column.

set -euo pipefail

in_file="$1"
out_file="$2"

awk -F'\t' '
NR==1 {print "CHR\tPosition\tChr\tbp\tA1\tA2\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8; next}
{
  split($1, a, ":")
  split($2, b, ":")
  print a[1] "\t" a[2] "\t" b[1] "\t" b[2] "\t" b[3] "\t" b[4] "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8
}' "$in_file" > "$out_file"
