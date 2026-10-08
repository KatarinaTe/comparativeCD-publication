#! /usr/bin/env bash

set -euo pipefail

# =============================================================================
# fetch-jars.sh
# Downloads the conform-gt and beagle jars used by this pipeline's CONFORM_GT
# and BEAGLE_IMPUTE processes, pinned to the exact releases run in the
# original 02_axiom_impute.sh.
# =============================================================================

BIN_DIR="bin"
mkdir -p "${BIN_DIR}"

wget -nc -P "${BIN_DIR}" "https://faculty.washington.edu/browning/conform-gt/conform-gt.24May16.cee.jar"
wget -nc -P "${BIN_DIR}" "https://faculty.washington.edu/browning/beagle/beagle.22Jul22.46e.jar"

chmod -R a-w "${BIN_DIR}"

echo "Jars ready in ${BIN_DIR}"
