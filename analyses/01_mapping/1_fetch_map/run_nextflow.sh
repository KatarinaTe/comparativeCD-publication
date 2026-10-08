#! /usr/bin/env bash

set -euo pipefail

trap 'echo "=!= ERROR: Workflow failed (exit code: $?) ==="' ERR

# Invoked once per batch by analyses/00_fetch-raw-data/run_batched.sh, which sets
# SAMPLESHEET (that batch's FASTQ samplesheet), WORK_DIR, and OUTPUT_DIR before calling
# this script. Sensible defaults are provided so this script also runs standalone
# (e.g. for the stub test) without the driver.
SAMPLESHEET="${SAMPLESHEET:-}"
WORK_DIR="${WORK_DIR:-work}"
OUTPUT_DIR="${OUTPUT_DIR:-../../../data/results/01_mapping/1_fetch_map/default}"

if [ -z "${SAMPLESHEET}" ]; then
    echo "=!= ERROR: SAMPLESHEET must be set (a batch's FASTQ samplesheet CSV) ===" >&2
    exit 1
fi

echo "=== Starting workflow: ==="
echo "=== Profile: ${NXF_PROFILE:-uppmax} ==="
echo "=== Samplesheet: ${SAMPLESHEET} ==="
WORKFLOW=../../../code/01_mapping/1_fetch_map/1_fetch_map.nf
nextflow run "${WORKFLOW}" \
    -params-file params.yml \
    --samplesheet "${SAMPLESHEET}" \
    -work-dir "${WORK_DIR}" \
    -output-dir "${OUTPUT_DIR}" \
    -profile "${NXF_PROFILE:-uppmax}" \
    ${CLI_OPTS:-""}

echo ""
echo "=== Cleaning up obsolete work directories ==="
nextflow clean -f -before last && find "${WORK_DIR}" -type d -empty -delete
echo "=== Workflow complete ==="
