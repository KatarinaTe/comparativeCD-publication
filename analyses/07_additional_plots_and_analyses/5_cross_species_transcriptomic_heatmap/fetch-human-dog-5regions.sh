#! /usr/bin/env bash

set -euo pipefail

# =============================================================================
# fetch-human-dog-5regions.sh
# Downloads human_dog_5regions.rds (the collaborator's integrated human-dog-mouse Seurat
# object, mouse cells already dropped; Figshare DOI 10.6084/m9.figshare.33787129) that
# cross_species_transcriptomic_heatmap.R consumes as external, already-computed input.
#
# Activated 2026-10-05: the deposit (same Figshare article as fetch-siletti-subset.sh's
# Siletti_human_subset_normalised.rds, now with this file added alongside it) went live.
# Confirmed via the Figshare API (api.figshare.com/v2/articles/33787129) and directly against
# the cluster-local copy fetched manually via gdown on 2026-09-30 (the one
# cross_species_transcriptomic_heatmap.R's development/verification run used, reproducing the
# published Extended Data Fig. 3 exactly) -- byte-identical, same size (5,115,368,146 bytes) and
# md5 (3812a1c0743bd506359fb77f1e158391). So the already-verified ED Fig. 3 run and this fetched
# copy are provably the same data, not just presumed to be.
# =============================================================================

BASE_DIR="../../../data/07_additional_plots_and_analyses"
mkdir -p "${BASE_DIR}"

# file_id:name:md5 — from the Figshare API (api.figshare.com/v2/articles/33787129)
FILE_ID="69486390"
NAME="human_dog_5regions.rds"
MD5="3812a1c0743bd506359fb77f1e158391"

dest="${BASE_DIR}/${NAME}"

if [ -f "${dest}" ] && [ "$(md5sum "${dest}" | cut -d' ' -f1)" = "${MD5}" ]; then
  echo "OK (already present, checksum verified): ${NAME}"
  exit 0
fi

wget -q -O "${dest}" "https://ndownloader.figshare.com/files/${FILE_ID}"

actual_md5="$(md5sum "${dest}" | cut -d' ' -f1)"
if [ "${actual_md5}" != "${MD5}" ]; then
  echo "CHECKSUM MISMATCH for ${NAME}: expected ${MD5}, got ${actual_md5}" >&2
  exit 1
fi
echo "Downloaded and verified: ${NAME}"
