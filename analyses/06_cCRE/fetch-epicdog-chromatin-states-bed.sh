#! /usr/bin/env bash

set -euo pipefail

# =============================================================================
# fetch-epicdog-chromatin-states-bed.sh
# Downloads the per-tissue EpicDog 13-state ChromHMM BED files that Fig. 4d's
# UU-vs-EpicDog comparison needs — already on CanFam4 coordinates, from the same
# Figshare record ("Dog Brain ChromHMM chromatin states files", DOI
# 10.6084/m9.figshare.31971747) that fetch-chromatin-states-bed.sh already pulls
# the UU/own 9-state files from.
#
# EpicDog's own publication (github.com/snu-cdrc/dog-reference-epigenome) called
# chromatin states on canFam3.1. This record's depositors did NOT liftOver that
# canFam3.1 output — they remapped EpicDog's raw ChIP/ATAC data to CanFam4 from
# scratch and re-ran peak calling + ChromHMM using the same ENCODE pipeline as the
# original. So there is no liftOver step here (an earlier version of this script
# fetched the canFam3.1 GitHub source and fed it through a canFam3.1->CanFam4
# liftOver pipeline at code/06_cCRE/liftover_epicdog/ — that was wrong and has been
# removed; see data/Manuscript_draft/REVIEW_AND_SUGGESTIONS.md:103).
#
# 13-state legend (from the Figshare record's description):
#   1  Bivalent TSS/Enh           8  Flanking active TSS2
#   2  Repressed polycomb         9  Active TSS
#   3  Repressed                 10  Quiescent
#   4  Heterochromatin           11  Active, weak enhancer
#   5  ZNF genes + repeats       12  Active, strong enhancer
#   6  Weak TSS                  13  Active, poised enhancer
#   7  Flanking active TSS1
# =============================================================================

BASE_DIR="../../data/06_cCRE/epicdog_chromatin_states_bed_canfam4"
mkdir -p "${BASE_DIR}"

# file_id:name:md5 — from the Figshare API (api.figshare.com/v2/articles/31971747)
FILES=(
  "68547067:cerebellum_13_dense.bed:b814c8711193e59dc3cc7c39fb0d75ab"
  "68547025:cerebrum_13_dense.bed:94a15f74b19a559b964d50508a1be606"
  "68547043:colon_13_dense.bed:22cff60f84aebf7862f5473364f35f96"
  "68547070:kidney_13_dense.bed:e80e46cf8c0839320666e79c56c5d001"
  "68547049:liver_13_dense.bed:520f4a8387cfaa8dd9fadccc913c4d98"
  "68547034:lung_13_dense.bed:053b75d3acf8b0dbabff251453083caa"
  "68547073:mammary_gland_13_dense.bed:9ee1747695a29265f3479a74bdc2daed"
  "68547037:ovary_13_dense.bed:09e52be85f6e333c4c7a21b720653cc3"
  "68547052:pancreas_13_dense.bed:fd81d5279b22939e00cf37cea365c646"
  "68547028:spleen_13_dense.bed:84376350b94d93997c57cc347f8238bf"
  "68547031:stomach_13_dense.bed:b11946eb268ed3e05f8f9639fc53b36a"
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

echo "EpicDog chromatin-state BED files (CanFam4) ready in ${BASE_DIR}"
