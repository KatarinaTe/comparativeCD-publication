#! /usr/bin/env bash

set -euo pipefail

# =============================================================================
# fetch-raw-data.sh
# Downloads reference files required for analyses/01_mapping/params.yml
#
# sha256 values below were computed directly against the files as hosted at
# kiddlabshare.med.umich.edu (2026-09-16) - this is a lab share site, not a
# versioned archive, so these checksums are what pin the content, not the URL.
# =============================================================================

BASE_DIR="../../data/raw-data"
REF_DIR="${BASE_DIR}/references"

mkdir -p "${REF_DIR}"

# fetch_and_verify <url> <dest> <expected_sha256>
# Skips the download if dest already exists and matches; always verifies after
# downloading, failing loudly on a mismatch rather than leaving a silently
# corrupt/stale file in place (replaces the old plain `wget -nc`).
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
# Genome reference (UU_Cfam_GSD_1.0_ROSY)
# ------------------------------------------------------------------------------
REF_BASE="https://kiddlabshare.med.umich.edu/public-data/UU_Cfam_GSD_1.0-Y"

fetch_and_verify "${REF_BASE}/UU_Cfam_GSD_1.0_ROSY.fa.gz" \
  "${REF_DIR}/UU_Cfam_GSD_1.0_ROSY.fa.gz" \
  "b359f60e54e4f51e6db5827acd99ae188d02b8db0c990b62b2858b4faaddc5aa"
gunzip -c "${REF_DIR}/UU_Cfam_GSD_1.0_ROSY.fa.gz" > "${REF_DIR}/UU_Cfam_GSD_1.0_ROSY.fa"

fetch_and_verify "${REF_BASE}/UU_Cfam_GSD_1.0_ROSY.fa.fai" \
  "${REF_DIR}/UU_Cfam_GSD_1.0_ROSY.fa.fai" \
  "be99071737a3d8867c30cac240058e5510c762326151496d43d1bff47270b382"
fetch_and_verify "${REF_BASE}/UU_Cfam_GSD_1.0_ROSY.dict" \
  "${REF_DIR}/UU_Cfam_GSD_1.0_ROSY.dict" \
  "6b80ea42b184f5d87a987ea12c4191501e57ea0dd5b030fc85833aa4dfe8dd31"

# ------------------------------------------------------------------------------
# BQSR known variants
# ------------------------------------------------------------------------------
fetch_and_verify "${REF_BASE}/UU_Cfam_GSD_1.0.BQSR.DB.bed.gz" \
  "${REF_DIR}/UU_Cfam_GSD_1.0.BQSR.DB.bed.gz" \
  "e01364f6ac10462eba3c32d3b40efec3425485c1f7491265b203648f4aafb95b"
fetch_and_verify "${REF_BASE}/UU_Cfam_GSD_1.0.BQSR.DB.bed.gz.tbi" \
  "${REF_DIR}/UU_Cfam_GSD_1.0.BQSR.DB.bed.gz.tbi" \
  "3987c085440ec7416568209291ef7c1c3ea05de790900a926479f346878142c9"

# ------------------------------------------------------------------------------
# Depth sites (CanineHD array)
# ------------------------------------------------------------------------------
fetch_and_verify "${REF_BASE}/SRZ189891_722g.simp.header.CanineHD.names.GSD_1.0.filter.vcf.gz" \
  "${REF_DIR}/SRZ189891_722g.simp.header.CanineHD.names.GSD_1.0.filter.vcf.gz" \
  "596797475123d2a510e22b91799f753a21be9a8059e7e01d7c0547de15d4a8f0"
fetch_and_verify "${REF_BASE}/SRZ189891_722g.simp.header.CanineHD.names.GSD_1.0.filter.vcf.gz.tbi" \
  "${REF_DIR}/SRZ189891_722g.simp.header.CanineHD.names.GSD_1.0.filter.vcf.gz.tbi" \
  "9da438d89f7720155b44d4bd3c4ba11100118d8aad57ab2f5e2e89d2d721778d"

# ------------------------------------------------------------------------------
# GATK imputation panel (phased-imputation-panel)
# ------------------------------------------------------------------------------
IMPUTE_BASE="https://kiddlabshare.med.umich.edu/dog10K/phased-imputation-panel"

fetch_and_verify "${IMPUTE_BASE}/AutoAndXPAR.Dog10K.phased.bcf" \
  "${REF_DIR}/AutoAndXPAR.Dog10K.phased.bcf" \
  "635e59f0f01d95fc87fc465f6b855e4c3d153b78eafb54416b64465844673927"
fetch_and_verify "${IMPUTE_BASE}/AutoAndXPAR.Dog10K.phased.bcf.csi" \
  "${REF_DIR}/AutoAndXPAR.Dog10K.phased.bcf.csi" \
  "ef343c174b8d7fe5b7fac222484cc1e7fd11f426bb0ab3189e984e437776e79b"

# ------------------------------------------------------------------------------
# Dog10K curated NCBI v106 gene models, and the gene-range file derived from them
# for `plink --clump-range` (needed by 04_gwas/1_mlma-loco and 2_polmm)
#
# genes6_UU_Cfam_GSD_1.0_ROSY.txt = every "gene" line of the GTF on chr1-38/X, as
# "chr start end gene_id". The awk below splits on whitespace like the original
# (2024-10-23), so the 278 pseudogenes whose source column is the two words
# "Curated Genomic" fall out exactly as they did then; verified byte-identical
# (sha256 below) to the file used for the published clumping.
# ------------------------------------------------------------------------------
ANNOT_BASE="https://kiddlabshare.med.umich.edu/dog10K/annotation"

fetch_and_verify "${ANNOT_BASE}/UU_Cfam_GSD_1.0_ROSY.refSeq.ensformat.gtf" \
  "${REF_DIR}/UU_Cfam_GSD_1.0_ROSY.refSeq.ensformat.gtf" \
  "ea8f86787cee42ee064cbc08eed18eca2e3997bd07add48b8be3f771bfeaaeb2"

GENE_RANGES="${BASE_DIR}/genes6_UU_Cfam_GSD_1.0_ROSY.txt"
awk '$3=="gene" {gsub(/"|;/,"",$10); sub(/^chr/,"",$1); if ($1 ~ /^([0-9]+|X)$/) print $1, $4, $5, $10}' \
  "${REF_DIR}/UU_Cfam_GSD_1.0_ROSY.refSeq.ensformat.gtf" > "${GENE_RANGES}"
actual="$(sha256sum "${GENE_RANGES}" | cut -d' ' -f1)"
if [ "${actual}" != "f6fa6ad1354ad5a5d056b8302ee83095fe86230f3fd0e0a85291d5de8f34bb2f" ]; then
  echo "CHECKSUM MISMATCH for $(basename "${GENE_RANGES}"): got ${actual}" >&2
  exit 1
fi
echo "Derived and verified: $(basename "${GENE_RANGES}")"

# ------------------------------------------------------------------------------
# canFam3 -> canFam4 liftover chain (needed by 02_data_processing/2_Axiom_imputation)
# ------------------------------------------------------------------------------
wget -nc -P "${REF_DIR}" "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/liftOver/canFam3ToCanFam4.over.chain.gz"

# ------------------------------------------------------------------------------
# Make the fetched reference data read-only (protects it from accidental
# edits/deletes) -- scoped to REF_DIR, not the whole of BASE_DIR: run_batched.sh
# later needs to create/write BASE_DIR/mapped/ (the permanent materialized-output
# archive) as a sibling of references/, and a BASE_DIR-wide lockdown blocks that
# mkdir for every user who runs the pipeline in its normal order.
# ------------------------------------------------------------------------------
chmod -R a-w "${REF_DIR}"

echo "Reference files ready in ${REF_DIR}"
