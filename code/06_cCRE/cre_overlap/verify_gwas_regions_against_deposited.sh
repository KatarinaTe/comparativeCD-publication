#!/bin/bash
set -euo pipefail

# verify_gwas_regions_against_deposited.sh -- one-off verification tool, NOT part of the live
# Nextflow DAG. Run this once 06b-clump-regions-100kb has been run for real (against real
# .clumped/.loco.mlma data) to check whether its live output actually matches the deposited
# reference copies cre_overlap/ currently uses in their place.
#
# Background: cre_overlap.nf's dogcd_gwas_regions/size_gwas_regions params point at two
# committed files (data/06_cCRE/dogCD_gws100kb_collapsed.bed,
# SIZE_gwas_clump100kb_regions_for_cCREintersect_clean.txt) that were designed to match what
# BUILD_CLUMP_REGIONS_100KB (06b) produces for the same inputs, but this has never actually been
# checked -- this environment has no real .clumped/.loco.mlma data to run 06b against (see
# clump_regions_100kb/process_definitions.nf's own header, and cre_overlap/README.md's Notes).
#
# Usage:
#   verify_gwas_regions_against_deposited.sh <06b_regions_dir>
# where <06b_regions_dir> is 06b's published output directory, containing
# dogCD_gws100kb_regions.bed and SIZE_gws100kb_regions.bed (see
# data/results/06_cCRE/clump_regions_100kb/regions/ after a real 06b run).
#
# What "pass" looks like: both diffs report zero added/removed base pairs. A pass means
# cre_overlap.nf's params.yml can be safely repointed at 06b's live output paths (and a real
# `depends-on = ["06b-clump-regions-100kb"]` added to the 06c-cre-overlap pixi task -- it's
# deliberately NOT there yet, see pixi.toml) instead of the deposited stand-ins used today.
# A fail means either the deposited copies are stale/wrong, or 06b's own build logic needs a
# closer look -- in either case, don't just switch the wiring without understanding why first.
#
# Requires: bedtools.

REGIONS_DIR="${1:?Usage: $0 <06b_regions_dir>}"

DEPOSITED_DOGCD="data/06_cCRE/dogCD_gws100kb_collapsed.bed"
DEPOSITED_SIZE="data/06_cCRE/SIZE_gwas_clump100kb_regions_for_cCREintersect_clean.txt"
LIVE_DOGCD="${REGIONS_DIR}/dogCD_gws100kb_regions.bed"
LIVE_SIZE="${REGIONS_DIR}/SIZE_gws100kb_regions.bed"

for f in "${DEPOSITED_DOGCD}" "${DEPOSITED_SIZE}" "${LIVE_DOGCD}" "${LIVE_SIZE}"; do
  [ -f "${f}" ] || { echo "=!= Missing: ${f}" >&2; exit 1; }
done

OUTDIR="$(mktemp -d)"

# Strip any leading '#'-comment header line (the deposited dogCD file has one; 06b's own
# bedtools sort|merge output never does), then sort+merge both sides identically so a byte
# comparison isn't thrown off by row order/formatting differences alone.
normalize() {
  grep -v '^#' "$1" | bedtools sort -i - | bedtools merge -i -
}

compare_one() {
  local label="$1" deposited="$2" live="$3"
  normalize "${deposited}" > "${OUTDIR}/${label}_deposited.bed"
  normalize "${live}"      > "${OUTDIR}/${label}_live.bed"

  local deposited_bp live_bp added_bp removed_bp
  deposited_bp=$(awk '{sum+=$3-$2} END{print sum+0}' "${OUTDIR}/${label}_deposited.bed")
  live_bp=$(awk '{sum+=$3-$2} END{print sum+0}' "${OUTDIR}/${label}_live.bed")
  # bp present in live but not deposited ("added"), and vice versa ("removed")
  added_bp=$(bedtools subtract -a "${OUTDIR}/${label}_live.bed" -b "${OUTDIR}/${label}_deposited.bed" \
    | awk '{sum+=$3-$2} END{print sum+0}')
  removed_bp=$(bedtools subtract -a "${OUTDIR}/${label}_deposited.bed" -b "${OUTDIR}/${label}_live.bed" \
    | awk '{sum+=$3-$2} END{print sum+0}')

  echo "--- ${label} ---"
  printf "%-20s %15s %15s\n" "" "deposited" "live (06b)"
  printf "%-20s %15s %15s\n" "total_bp" "${deposited_bp}" "${live_bp}"
  printf "%-20s %15s\n" "bp only in live"      "${added_bp}"
  printf "%-20s %15s\n" "bp only in deposited" "${removed_bp}"
  if [ "${added_bp}" -eq 0 ] && [ "${removed_bp}" -eq 0 ]; then
    echo "PASS: identical region sets"
  else
    echo "FAIL: region sets differ -- do not switch cre_overlap.nf's wiring without understanding why"
  fi
  echo
}

compare_one "dogCD" "${DEPOSITED_DOGCD}" "${LIVE_DOGCD}"
compare_one "SIZE"  "${DEPOSITED_SIZE}"  "${LIVE_SIZE}"
