#! /usr/bin/env bash

set -euo pipefail

# =============================================================================
# fetch-siletti-subset.sh
# Downloads the Siletti human brain snRNA-seq subset (Figshare DOI
# 10.6084/m9.figshare.33787129) that cross_species_expression.R consumes as external,
# already-computed input (the collaborator confirmed the script that produced this subset
# from the full Siletti et al. 2023 atlas cannot be provided — see REPRODUCIBILITY_AUDIT.md).
# =============================================================================

BASE_DIR="../../../data/07_additional_plots_and_analyses"
mkdir -p "${BASE_DIR}"

# file_id:name:md5 — from the Figshare API (api.figshare.com/v2/articles/33787129)
FILE_ID="68624272"
NAME="Siletti_human_subset_normalised.rds"
MD5="9386b23e001dd0e09c9d5b89b9935fe6"

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
