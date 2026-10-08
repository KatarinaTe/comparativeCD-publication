#! /usr/bin/env bash

set -euo pipefail

# =============================================================================
# fetch-scilifelab-deposit.sh
# Downloads the pipeline inputs deposited on the SciLifeLab Data Repository,
# https://doi.org/10.17044/scilifelab.33339309 (version 4, CC BY 4.0), into the
# paths the pipeline's params.yml files expect.
#
# Files are fetched by their Figshare file IDs (stable within a deposit version).
# sha256 values were computed against the files as downloaded on 2026-10-08; each
# file's md5 also matches the md5 the repository itself reports for it.
#
# Not fetched here (not pipeline inputs; the pipeline regenerates them): the merged
# DA_MERGED_GENCOVE_AXIOM_QC5 plinkset (9.7 GB), the per-trait GWAS summary
# statistics and ALLFAM_sex_age_phenos_batch.txt. Download those from the DOI above
# if you want the published versions.
# =============================================================================

FIGSHARE="https://ndownloader.figshare.com/files"
AXIOM_DIR="../../data/raw-data/axiom"
DP_DIR="../../data/02_data_processing"

mkdir -p "${AXIOM_DIR}" "${DP_DIR}"

# fetch_and_verify <url> <dest> <expected_sha256> -- same behaviour as in
# fetch-reference-data.sh: skip if already present and matching, fail loudly on a
# mismatch.
fetch_and_verify() {
  local url="$1" dest="$2" expected="$3" actual

  if [ -f "${dest}" ]; then
    actual="$(sha256sum "${dest}" | cut -d' ' -f1)"
    if [ "${actual}" = "${expected}" ]; then
      echo "OK (already present, checksum verified): $(basename "${dest}")"
      return
    fi
  fi

  wget -q -O "${dest}" "${url}"

  actual="$(sha256sum "${dest}" | cut -d' ' -f1)"
  if [ "${actual}" != "${expected}" ]; then
    echo "CHECKSUM MISMATCH for $(basename "${dest}"): expected ${expected}, got ${actual}" >&2
    exit 1
  fi
  echo "Downloaded and verified: $(basename "${dest}")"
}

# ------------------------------------------------------------------------------
# Axiom Canine Genotyping Array A+B genotypes, canFam3 (02b-axiom-imputation input)
#
# Deposited as affy_round3.{bed,bim,fam}: 411 dogs x 1,011,992 variants. This is
# the original 804-dog Affy_merged plinkset already restricted to the 411 dogs in
# data/02_data_processing/ind_to_keep_round3.txt (identical sample IDs) -- i.e. the
# output of LIFTOVER_TO_CANFAM4's own first `plink --keep` step, so re-applying that
# filter is a no-op. Saved under a different prefix because that step writes its
# own output as affy_round3.*.
# ------------------------------------------------------------------------------
fetch_and_verify "${FIGSHARE}/67885554" "${AXIOM_DIR}/axiom_411_canfam3.bed" \
  "d08b89bb039ba92cffce440d4dca2263f6522e655fda057919b48b0f7e7083c9"
fetch_and_verify "${FIGSHARE}/67885551" "${AXIOM_DIR}/axiom_411_canfam3.bim" \
  "77e830bcc10f86c49a96b461567a3b1f78b678989455e1081a356e19f280b7e6"
fetch_and_verify "${FIGSHARE}/67885548" "${AXIOM_DIR}/axiom_411_canfam3.fam" \
  "ad912dcc02ae8d06a85cc7655068b4021a49a50749a922aa4718bfa05bb58afd"

# ------------------------------------------------------------------------------
# Darwin's Ark survey and metadata subsets (02c-survey-data, 02e-plotting inputs)
# ------------------------------------------------------------------------------
fetch_and_verify "${FIGSHARE}/69717552" "${DP_DIR}/response_df_dogCDitems.txt" \
  "bec2b6da56b7e8835efa85acf840581c32dde5c5ef3e4001a03e00ee3cb1edc0"
fetch_and_verify "${FIGSHARE}/69678312" "${DP_DIR}/DarwinsArk_20220715_dogs_genotyped_breed_sex.csv" \
  "1356fd4bab0597b3dfb7a494fd0202f2c7eb91aabb2d3bb3273c74135481bebc"

echo "SciLifeLab deposit files ready (Axiom: ${AXIOM_DIR}; survey/metadata: ${DP_DIR})"
