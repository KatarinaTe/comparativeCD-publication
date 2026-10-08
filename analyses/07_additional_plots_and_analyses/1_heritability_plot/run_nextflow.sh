#! /usr/bin/env bash

set -euo pipefail

trap 'echo "=!= ERROR: Workflow failed (exit code: $?) ==="' ERR

echo "=== Starting workflow: ==="
echo "=== Profile: ${NXF_PROFILE:-uppmax} ==="
WORKFLOW=../../../code/07_additional_plots_and_analyses/1_heritability_plot/1_heritability_plot.nf
nextflow run "${WORKFLOW}" \
    -params-file params.yml \
    -profile "${NXF_PROFILE:-uppmax}" \
    ${CLI_OPTS:-""}

echo ""
echo "=== Cleaning up obsolete work directories ==="
nextflow clean -f -before last && find work -type d -empty -delete
echo "=== Workflow complete ==="
