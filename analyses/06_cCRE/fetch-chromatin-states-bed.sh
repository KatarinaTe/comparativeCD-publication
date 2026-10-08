#! /usr/bin/env bash

set -euo pipefail

# =============================================================================
# fetch-chromatin-states-bed.sh
# Downloads the per-brain-region ChromHMM 9-state BED files ("Dog Brain ChromHMM
# chromatin state bed files", Figshare DOI 10.6084/m9.figshare.31971747) that
# code/06_cCRE and code/04_gwas/5_gwas_regions consume as external, already-computed
# input (they are the state.bed outputs behind Fig. 4c/d of the manuscript).
# =============================================================================

BASE_DIR="../../data/06_cCRE/chromatin_states_bed"
mkdir -p "${BASE_DIR}"

# file_id:name:md5 — from the Figshare API (api.figshare.com/v2/articles/31971747)
FILES=(
  "63584772:ACG_9_dense.bed:8f839e56e0916d6f400cfb5ca86cd954"
  "63584775:cerebellum_9_dense.bed:068aed094493091e42dd55b3b2a373e0"
  "63584781:frontal_lobe_9_dense.bed:ebb633298a0f421d4b01566c195a9d57"
  "63584778:hypothalamus_9_dense.bed:c28fc981e567506d54cf5c4ef5c69bc2"
  "63584787:occipital_cortex_9_dense.bed:b3d4699823118a7b9ae97143e201ec38"
  "63584784:striatum_9_dense.bed:b6f6b8a0b22f653fd45408ddcf78b642"
  "63584790:temporal_cortex_9_dense.bed:f7aecda8eef8ee9d595951c225948426"
  "63584793:thalamus_9_dense.bed:2a685d24ab9ddb539ca01c1f021779ee"
)

for entry in "${FILES[@]}"; do
  IFS=":" read -r file_id name md5 <<< "${entry}"
  dest="${BASE_DIR}/${name}"

  if [ -f "${dest}" ] && [ "$(md5sum "${dest}" | cut -d' ' -f1)" = "${md5}" ]; then
    echo "OK (already present, checksum verified): ${name}"
    continue
  fi

  wget -q -O "${dest}" "https://ndownloader.figshare.com/files/${file_id}"

  actual_md5="$(md5sum "${dest}" | cut -d' ' -f1)"
  if [ "${actual_md5}" != "${md5}" ]; then
    echo "CHECKSUM MISMATCH for ${name}: expected ${md5}, got ${actual_md5}" >&2
    exit 1
  fi
  echo "Downloaded and verified: ${name}"
done

echo "Chromatin-state BED files ready in ${BASE_DIR}"
