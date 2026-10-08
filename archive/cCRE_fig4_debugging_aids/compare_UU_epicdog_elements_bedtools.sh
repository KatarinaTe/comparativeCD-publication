#!/bin/bash
set -euo pipefail

# Cross-check for compare_UU_epicdog_elements.R (Fig. 4d): reimplements the same logic
# in pure bedtools, so the R (valr) and bedtools results can be compared directly to
# check whether a discrepancy comes from a valr-vs-bedtools implementation difference
# or from something else (e.g. how all_tissues_filtered.bed was actually built).
#
# Requires: bedtools (env/06_cCRE/clump_regions.yml's container, or `module load bedtools`).
#
# Inputs (fetch first if not already present):
#   analyses/06_cCRE/fetch-chromatin-states-bed.sh          -> data/06_cCRE/chromatin_states_bed/
#   analyses/06_cCRE/fetch-epicdog-chromatin-states-bed.sh   -> data/06_cCRE/epicdog_chromatin_states_bed_canfam4/

UU_DIR="${UU_DIR:-data/06_cCRE/chromatin_states_bed}"
EPIC_DIR="${EPIC_DIR:-data/06_cCRE/epicdog_chromatin_states_bed_canfam4}"
OUTDIR="${OUTDIR:-./tmp_fig4d_bedtools_check}"
mkdir -p "${OUTDIR}"

strip_track() {
  # Drop the UCSC "track name=..." header line, same fix as read_dense_bed() in R
  tail -n +2 "$1"
}

# --- Step 1: per-state merge across all 8 UU regions (state 7 simply excluded) ------
echo "Step 1: building all_tissues_filtered.bed (per-state merge, no combined states)..."
> "${OUTDIR}/all_tissues_filtered.bed"
for state in 1 2 3 4 5 6 8 9; do
  for f in "${UU_DIR}"/*_9_dense.bed; do
    strip_track "$f"
  done \
    | awk -v s="$state" 'BEGIN{FS="\t";OFS="\t"} $4==s {print $1,$2,$3}' \
    | bedtools sort -i - \
    | bedtools merge -i - \
    | awk -v s="$state" 'BEGIN{OFS="\t"} {print $1,$2,$3,s}' \
    >> "${OUTDIR}/all_tissues_filtered.bed"
done
echo "  $(wc -l < "${OUTDIR}/all_tissues_filtered.bed") rows"

# --- Step 2: EPIC cerebellum+cerebrum, header-stripped, ready to filter/merge below -
strip_track "${EPIC_DIR}/cerebellum_13_dense.bed" > "${OUTDIR}/epic_cerebellum.bed"
strip_track "${EPIC_DIR}/cerebrum_13_dense.bed"   > "${OUTDIR}/epic_cerebrum.bed"

# --- Step 3: per-category comparison, any-overlap boolean classification (matches
#     the R script's semi/anti-join logic: -u = "has >=1 overlap" (Shared),
#     -v = "has 0 overlap" (Unique)) ---------------------------------------------
compare_category() {
  local label="$1" epic_states="$2" uu_states="$3"

  awk -v states="$epic_states" 'BEGIN{FS="\t";OFS="\t"; n=split(states,a,","); for(i=1;i<=n;i++) want[a[i]]=1}
       want[$4] {print $1,$2,$3}' "${OUTDIR}/epic_cerebellum.bed" "${OUTDIR}/epic_cerebrum.bed" \
    | bedtools sort -i - | bedtools merge -i - > "${OUTDIR}/epic_${label}.bed"

  awk -v states="$uu_states" 'BEGIN{FS="\t";OFS="\t"; n=split(states,a,","); for(i=1;i<=n;i++) want[a[i]]=1}
       want[$4] {print $1,$2,$3}' "${OUTDIR}/all_tissues_filtered.bed" \
    | bedtools sort -i - | bedtools merge -i - > "${OUTDIR}/uu_${label}.bed"

  n_epic=$(wc -l < "${OUTDIR}/epic_${label}.bed")
  n_uu=$(wc -l < "${OUTDIR}/uu_${label}.bed")
  shared=$(bedtools intersect -a "${OUTDIR}/epic_${label}.bed" -b "${OUTDIR}/uu_${label}.bed" -u | wc -l)
  unique_epic=$(bedtools intersect -a "${OUTDIR}/epic_${label}.bed" -b "${OUTDIR}/uu_${label}.bed" -v | wc -l)
  unique_uu=$(bedtools intersect -a "${OUTDIR}/uu_${label}.bed" -b "${OUTDIR}/epic_${label}.bed" -v | wc -l)

  echo -e "${label}\t${unique_epic}\t${unique_uu}\t${shared}\t(epic total=${n_epic}, uu total=${n_uu})"
}

echo ""
echo "Feature	Unique_Epic	Unique_UU	Shared"
compare_category "Promoters" "6,7,8,9" "1,4,5,6"
compare_category "Enhancers" "11,12,13" "2,8"
compare_category "Repressed" "1,2,3" "9"
