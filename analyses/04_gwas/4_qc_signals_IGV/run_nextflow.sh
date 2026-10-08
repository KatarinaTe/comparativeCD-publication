#! /usr/bin/env bash

set -euo pipefail

trap 'echo "=!= ERROR: Workflow failed (exit code: $?) ==="' ERR

echo "=== Starting workflow: ==="
echo "=== Profile: ${NXF_PROFILE:-uppmax} ==="
WORKFLOW=../../../code/04_gwas/4_qc_signals_IGV/4_qc_signals_IGV.nf
nextflow run "${WORKFLOW}" \
    -params-file params.yml \
    -profile "${NXF_PROFILE:-uppmax}" \
    ${CLI_OPTS:-""}

echo ""
echo "=== Cleaning up obsolete work directories ==="
nextflow clean -f -before last && find work -type d -empty -delete
echo "=== Workflow complete ==="
