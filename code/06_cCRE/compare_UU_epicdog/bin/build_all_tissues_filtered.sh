#!/usr/bin/env bash
set -euo pipefail

# Fig. 4d's UU-side input for compare_UU_epicdog_elements.R: concatenates the 8
# per-brain-region UU 9-state ChromHMM BEDs, tags each row with its source region, and
# drops unannotated (state 7) rows.
#
# Wired into Nextflow 2026-09-23 -- reads/writes directly in the task work dir, since
# process_definitions.nf stages each input under its own real basename (ACG_9_dense.bed,
# etc.). The loose original (code/06_cCRE/build_all_tissues_filtered.sh, SRC_DIR/OUTPUT env
# vars relative to analyses/06_cCRE) is left untouched for standalone/manual use.
#
# Provenance, confirmed directly by the manuscript author (2026-09-21):
# "The 'all_tissues_filtered.bed' file is a concatenated file of all the
# genomic regions across all the brain regions we sampled... I simply
# concatenated all the bed files for each brain region into a single bed
# file, and then filtered out all regions labelled state 7 (unannotated),
# so only annotated states are represented. Column 4 containing the states
# should only ever contain one value as each row is a specific state
# identified in a specific tissue."

OUTPUT="all_tissues_filtered.bed"

# filename:region-label pairs -- region labels match the real deposited file's own
# column-5 spelling exactly (e.g. "frontal_cortex", not "frontal_lobe", even though the
# source filename says frontal_lobe).
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

> "${OUTPUT}"
for pair in "${PAIRS[@]}"; do
  IFS=":" read -r f region <<< "${pair}"
  # NR>1: drop this file's own UCSC track header line. $4!=7: drop unannotated-state rows.
  # Tag with this file's brain region as column 5.
  awk -F'\t' -v OFS='\t' -v region="${region}" \
    'NR>1 && $4 != 7 {print $1,$2,$3,$4,region}' "${f}" >> "${OUTPUT}"
done

echo "Wrote $(wc -l < "${OUTPUT}") rows to ${OUTPUT}"
