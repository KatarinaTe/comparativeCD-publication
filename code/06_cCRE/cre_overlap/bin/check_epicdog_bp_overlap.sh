#!/bin/bash
set -euo pipefail

# Preliminary test script for the EpicDog tissue-specificity analysis (formerly Fig. 4e; "Odds ratios for overlap of dogCD and
# dogSize GWAS regions with brain cCREs versus non-brain tissues (EpicDog
# resource)"). Computes the bp-level overlap numbers that
# testing_CRE_overlap_260617.R's `tissue_specificity_plot` figure currently
# hardcodes (dogCD_O_brain, dogCD_O_other, SIZE_O_brain, SIZE_O_other,
# brain_cre_bp, other_cre_bp), so they can be recomputed live instead.
#
# Combines check_CRE_overlap.sh's brain-vs-other tissue pooling logic with
# SIZE_check_UU_CRE_overlap.sh's bp-sum-after-intersect logic (both already in
# this directory), applied to the EpicDog per-tissue BEDs instead of the UU
# per-region ones.
#
# Requires: bedtools (already frozen in env/06_cCRE/clump_regions.yml's
# container, community.wave.seqera.io/library/clump_regions:09876be1daa6bcb4 —
# or `module load bedtools` on the HPC).
#
# Inputs (run analyses/06_cCRE/fetch-epicdog-chromatin-states-bed.sh first if
# not already present):
EPICDOG_DIR="${EPICDOG_DIR:-data/06_cCRE/epicdog_chromatin_states_bed_canfam4}"
DOGCD_GWAS="${DOGCD_GWAS:-data/06_cCRE/dogCD_gws100kb_collapsed.bed}"
SIZE_GWAS="${SIZE_GWAS:-data/06_cCRE/SIZE_gwas_clump100kb_regions_for_cCREintersect_clean.txt}"
OUTDIR="${OUTDIR:-$(mktemp -d)}"
OUTPUT="${OUTPUT:-data/06_cCRE/epicdog_tissue_bp_overlaps.txt}"

BRAIN_TISSUES=(cerebellum cerebrum)
OTHER_TISSUES=(colon kidney liver lung mammary_gland ovary pancreas spleen stomach)

echo "Working directory for intermediate files: ${OUTDIR}"
mkdir -p "${OUTDIR}"

# Step 1: filter each tissue's 13-state BED to promoter-like (7,8,9) and
# enhancer-like (11,12,13) states — same states check_CRE_overlap.sh's own
# element-count version already uses.
echo "Step 1: filtering active states (7,8,9,11,12,13) per tissue..."
for tissue in "${BRAIN_TISSUES[@]}" "${OTHER_TISSUES[@]}"; do
  awk 'BEGIN{FS="\t";OFS="\t"} $4 ~ /^(7|8|9|11|12|13)$/ {print $1,$2,$3,$4}' \
    "${EPICDOG_DIR}/${tissue}_13_dense.bed" > "${OUTDIR}/${tissue}_active.bed"
done

# Step 2: pool brain (cerebellum+cerebrum) vs. other (9 tissues), sort+merge
# to a set of non-overlapping active-element regions per category.
echo "Step 2: pooling and merging brain vs. other tissue categories..."
cat "${OUTDIR}"/cerebellum_active.bed "${OUTDIR}"/cerebrum_active.bed | \
  bedtools sort -i - | bedtools merge -i - > "${OUTDIR}/brain_merged.bed"

cat $(for t in "${OTHER_TISSUES[@]}"; do echo "${OUTDIR}/${t}_active.bed"; done) | \
  bedtools sort -i - | bedtools merge -i - > "${OUTDIR}/other_merged.bed"

brain_cre_bp=$(awk '{sum+=$3-$2} END{print sum+0}' "${OUTDIR}/brain_merged.bed")
other_cre_bp=$(awk '{sum+=$3-$2} END{print sum+0}' "${OUTDIR}/other_merged.bed")

# Step 3: strip the dogCD file's header comment (bedtools should skip '#'
# lines automatically, but do it explicitly to be safe/deterministic).
grep -v '^#' "${DOGCD_GWAS}" > "${OUTDIR}/dogCD_gwas.bed"
cp "${SIZE_GWAS}" "${OUTDIR}/SIZE_gwas.bed"

# Step 4: intersect each GWAS region set against brain_merged/other_merged,
# sort+merge the overlap segments (avoids double-counting if GWAS regions
# overlap each other), then sum bp — same approach as
# SIZE_check_UU_CRE_overlap.sh.
bp_overlap() {
  local elements_bed="$1" gwas_bed="$2"
  bedtools intersect -a "${elements_bed}" -b "${gwas_bed}" | \
    bedtools sort -i - | bedtools merge -i - | \
    awk '{sum+=$3-$2} END{print sum+0}'
}

echo "Step 3: intersecting with dogCD and SIZE GWAS regions..."
dogCD_O_brain=$(bp_overlap "${OUTDIR}/brain_merged.bed" "${OUTDIR}/dogCD_gwas.bed")
dogCD_O_other=$(bp_overlap "${OUTDIR}/other_merged.bed" "${OUTDIR}/dogCD_gwas.bed")
SIZE_O_brain=$(bp_overlap "${OUTDIR}/brain_merged.bed" "${OUTDIR}/SIZE_gwas.bed")
SIZE_O_other=$(bp_overlap "${OUTDIR}/other_merged.bed" "${OUTDIR}/SIZE_gwas.bed")

echo "--------------------------------------------------"
echo "Computed vs. currently-hardcoded values (testing_CRE_overlap_260617.R lines 358-364):"
echo "--------------------------------------------------"
printf "%-16s %15s %15s\n" "" "computed" "hardcoded"
printf "%-16s %15s %15s\n" "brain_cre_bp"  "${brain_cre_bp}"  "247096600"
printf "%-16s %15s %15s\n" "other_cre_bp"  "${other_cre_bp}"  "393644600"
printf "%-16s %15s %15s\n" "dogCD_O_brain" "${dogCD_O_brain}" "799607"
printf "%-16s %15s %15s\n" "dogCD_O_other" "${dogCD_O_other}" "1109740"
printf "%-16s %15s %15s\n" "SIZE_O_brain"  "${SIZE_O_brain}"  "710315"
printf "%-16s %15s %15s\n" "SIZE_O_other"  "${SIZE_O_other}"  "1330316"
echo "--------------------------------------------------"

# Write a committed-style intermediate file, same shape as the UU files
# (data/06_cCRE/{brain,SIZE}_region_bp_overlaps.txt), so the R script can
# eventually read this instead of the hardcoded scalars.
{
  echo -e "GWAS_set\tTissue_category\tOverlap_bp\tTotal_bp"
  echo -e "dogCD\tbrain\t${dogCD_O_brain}\t${brain_cre_bp}"
  echo -e "dogCD\tother\t${dogCD_O_other}\t${other_cre_bp}"
  echo -e "SIZE\tbrain\t${SIZE_O_brain}\t${brain_cre_bp}"
  echo -e "SIZE\tother\t${SIZE_O_other}\t${other_cre_bp}"
} > "${OUTPUT}"
echo "Wrote ${OUTPUT}"
