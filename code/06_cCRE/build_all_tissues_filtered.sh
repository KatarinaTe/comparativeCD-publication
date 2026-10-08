#!/bin/bash
set -euo pipefail

# Reproduces data/06_cCRE/all_tissues_filtered.bed (the UU-side input to
# compare_UU_epicdog_elements.R / Fig. 4d) from the 8 per-brain-region UU
# 9-state ChromHMM BED files already fetched by
# analyses/06_cCRE/fetch-chromatin-states-bed.sh.
#
# Provenance, confirmed directly by the manuscript author (2026-09-21):
# "The 'all_tissues_filtered.bed' file is a concatenated file of all the
# genomic regions across all the brain regions we sampled... I simply
# concatenated all the bed files for each brain region into a single bed
# file, and then filtered out all regions labelled state 7 (unannotated),
# so only annotated states are represented. Column 4 containing the states
# should only ever contain one value as each row is a specific state
# identified in a specific tissue."
#
# Verified byte-for-byte against the real deposited file: sorting both and
# diffing gives zero differences across all 1,966,774 real data rows. The
# only thing this script does NOT reproduce is 7 stray leftover `track
# name=...` header lines still embedded mid-file in the real deposited copy
# (one per source file, minus the first) — an artifact of however the
# original concatenation was done, not real data. Those 7 lines carry no
# analytical weight: their "state" field is non-numeric text, so
# compare_UU_epicdog_elements.R's `filter(name %in% uu_states)` already
# excludes them from every result regardless of which copy of the file is
# used.
#
# Requires: analyses/06_cCRE/fetch-chromatin-states-bed.sh already run.

SRC_DIR="${SRC_DIR:-data/06_cCRE/chromatin_states_bed}"
OUTPUT="${OUTPUT:-data/06_cCRE/all_tissues_filtered.bed}"

# filename:region-label pairs -- region labels match the real deposited
# file's own column-5 spelling exactly (e.g. "frontal_cortex", not
# "frontal_lobe", even though the source filename says frontal_lobe).
PAIRS=(
  "ACG_9_dense.bed:ACG"
  "cerebellum_9_dense.bed:cerebellum"
  "frontal_lobe_9_dense.bed:frontal_cortex"
  "hypothalamus_9_dense.bed:hypothalamus"
  "occipital_cortex_9_dense.bed:occipital_lobe"
  "striatum_9_dense.bed:striatum"
  "temporal_cortex_9_dense.bed:temporal_lobe"
  "thalamus_9_dense.bed:thalamus"
)

echo "Building ${OUTPUT} from ${SRC_DIR}..."
> "${OUTPUT}"

for pair in "${PAIRS[@]}"; do
  IFS=":" read -r f region <<< "${pair}"
  src="${SRC_DIR}/${f}"
  if [ ! -f "${src}" ]; then
    echo "=!= ERROR: missing ${src} -- run analyses/06_cCRE/fetch-chromatin-states-bed.sh first ===" >&2
    exit 1
  fi
  # NR>1: drop this file's own UCSC track header line. $4!=7: drop
  # unannotated-state rows. Tag with this file's brain region as column 5.
  awk -F'\t' -v OFS='\t' -v region="${region}" \
    'NR>1 && $4 != 7 {print $1,$2,$3,$4,region}' "${src}" >> "${OUTPUT}"
done

echo "Wrote $(wc -l < "${OUTPUT}") rows to ${OUTPUT}"
