#! /usr/bin/env bash

set -euo pipefail

trap 'echo "=!= ERROR: Workflow failed (exit code: $?) ==="' ERR

# Run exactly once, after every 1_fetch_map batch (both tracks) has completed and been
# materialized — see analyses/00_fetch-raw-data/run_batched.sh. Running this more than
# once per cohort, or before every batch has landed, silently changes the joint-calling
# and GLIMPSE imputation results — see code/01_mapping/2_cohort_genotyping/2_cohort_genotyping.nf's
# header and archive/conversion_notes/01_mapping.md.

echo "=== Starting workflow: ==="
echo "=== Profile: ${NXF_PROFILE:-uppmax} ==="
WORKFLOW=../../../code/01_mapping/2_cohort_genotyping/2_cohort_genotyping.nf
nextflow run "${WORKFLOW}" \
    -params-file params.yml \
    -profile "${NXF_PROFILE:-uppmax}" \
    ${CLI_OPTS:-""}

echo ""
echo "=== Cleaning up obsolete work directories ==="
nextflow clean -f -before last && find work -type d -empty -delete
echo "=== Workflow complete ==="
