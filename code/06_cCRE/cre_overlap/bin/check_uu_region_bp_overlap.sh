#!/bin/bash
set -euo pipefail

# check_uu_region_bp_overlap.sh -- bp-level overlap of dogCD and SIZE GWAS regions with each of
# the 8 UU dog-brain-region chromatin-state BEDs, for Fig. 4e (formerly 4f) ("Odds ratios for dogCD versus
# dogSize GWAS region overlap with cCREs across eight dog brain regions (UU dataset)").
#
# Replaces two things that predate this script: SIZE_check_UU_CRE_overlap.sh (SIZE only,
# hardcoded absolute UPPMAX paths, and depended on per-region *_UU_CRE_active_elements.bed files
# from a commented-out, never-run generation loop) and a dogCD-side script that was never
# committed at all -- brain_region_bp_overlaps.txt (dogCD) has been a deposited-only result up to
# now. This single script computes both from the same source chromatin-state BEDs used elsewhere
# in this repo (data/06_cCRE/chromatin_states_bed/*_9_dense.bed), the same active-state definition
# SIZE_check_UU_CRE_overlap.sh's own (commented-out) filtering step used (promoter-like 1,4,5,6;
# enhancer-like 2,8; open chromatin 3 -- i.e. every state except 7 "no signal" and 9 "repressed").
#
# Verified 2026-09-21: recomputing this way reproduces the currently-deposited
# brain_region_bp_overlaps.txt and SIZE_brain_region_bp_overlaps.txt exactly, region-for-region,
# byte-for-byte -- no equivalent of the mammary_gland/mammary bug found here (this analysis has no
# tissue-pooling step to have that class of bug in the first place).
#
# Requires: bedtools (already frozen in env/06_cCRE/clump_regions.yml's container,
# community.wave.seqera.io/library/clump_regions:09876be1daa6bcb4 -- or `module load bedtools`).

STATES_DIR="${STATES_DIR:-data/06_cCRE/chromatin_states_bed}"
DOGCD_GWAS="${DOGCD_GWAS:-data/06_cCRE/dogCD_gws100kb_collapsed.bed}"
SIZE_GWAS="${SIZE_GWAS:-data/06_cCRE/SIZE_gwas_clump100kb_regions_for_cCREintersect_clean.txt}"
OUTDIR="${OUTDIR:-$(mktemp -d)}"
DOGCD_OUTPUT="${DOGCD_OUTPUT:-data/06_cCRE/brain_region_bp_overlaps.txt}"
SIZE_OUTPUT="${SIZE_OUTPUT:-data/06_cCRE/SIZE_brain_region_bp_overlaps.txt}"

echo "Working directory for intermediate files: ${OUTDIR}"
mkdir -p "${OUTDIR}"

grep -v '^#' "${DOGCD_GWAS}" > "${OUTDIR}/dogCD_gwas.bed"
cp "${SIZE_GWAS}" "${OUTDIR}/SIZE_gwas.bed"

# Region name -> source BED filename (9-state UU chromatin-state calls, one file per region).
# Parallel arrays, not an associative array, so this runs under plain bash 3.2 too (e.g. macOS).
REGIONS=(ACG cerebellum frontal hypothalamus occipital striatum temporal thalamus)
REGION_FILES=(ACG_9_dense.bed cerebellum_9_dense.bed frontal_lobe_9_dense.bed hypothalamus_9_dense.bed occipital_cortex_9_dense.bed striatum_9_dense.bed temporal_cortex_9_dense.bed thalamus_9_dense.bed)

echo -e "Region\tOverlap_bp\tTotal_bp" > "${DOGCD_OUTPUT}"
echo -e "Region\tOverlap_bp\tTotal_bp" > "${SIZE_OUTPUT}"

for i in "${!REGIONS[@]}"; do
  region="${REGIONS[$i]}"
  src="${STATES_DIR}/${REGION_FILES[$i]}"
  echo "Processing ${region}..."

  # Filter to active states (exclude 7 "no signal" and 9 "repressed"); skip the UCSC track header line.
  awk 'BEGIN{FS="\t";OFS="\t"} $4 ~ /^(1|2|3|4|5|6|8)$/ {print $1,$2,$3,$4}' "${src}" \
    > "${OUTDIR}/${region}_active.bed"

  total_bp=$(awk '{sum+=$3-$2} END{print sum+0}' "${OUTDIR}/${region}_active.bed")

  dogcd_overlap=$(bedtools intersect -a "${OUTDIR}/${region}_active.bed" -b "${OUTDIR}/dogCD_gwas.bed" | \
    bedtools sort -i - | bedtools merge -i - | awk '{sum+=$3-$2} END{print sum+0}')
  size_overlap=$(bedtools intersect -a "${OUTDIR}/${region}_active.bed" -b "${OUTDIR}/SIZE_gwas.bed" | \
    bedtools sort -i - | bedtools merge -i - | awk '{sum+=$3-$2} END{print sum+0}')

  echo -e "${region}\t${dogcd_overlap}\t${total_bp}" >> "${DOGCD_OUTPUT}"
  echo -e "${region}\t${size_overlap}\t${total_bp}" >> "${SIZE_OUTPUT}"
done

echo "Wrote ${DOGCD_OUTPUT} and ${SIZE_OUTPUT}"
