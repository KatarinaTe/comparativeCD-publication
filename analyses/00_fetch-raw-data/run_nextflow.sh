#! /usr/bin/env bash
#SBATCH -A uppmax2025-2-420
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --mem=16G
#SBATCH -t 2-00:00:00
#SBATCH -J fetch-raw-data
#SBATCH --chdir=.

# Runs as a full-node SLURM batch job on Pelle so the Nextflow head process
# (and its polling of per-task state) runs on a compute node instead of
# relying on many individually-tracked slurm/uppmax jobs. Submit with:
#   sbatch run_nextflow.sh

set -euo pipefail

trap 'echo "=!= ERROR: Workflow failed (exit code: $?) ==="' ERR

# Anchor to the submission directory: SLURM does not guarantee the batch
# job's initial working directory matches where it was submitted from, and
# the script bytes SLURM actually executes are a spooled copy, so "$0"/
# BASH_SOURCE can't be used to recover the original path. Fall back to the
# script's own directory for direct (non-sbatch) invocation.
cd "${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"

echo "=== Starting workflow: ==="
echo "=== Profile: ${NXF_PROFILE:-apptainer} ==="
export NXF_SYNTAX_PARSER=v1
WORKFLOW=nf-core/fetchngs
nextflow run "${WORKFLOW}" \
    -r 1.12.0 \
    -params-file params.yml \
    -profile "${NXF_PROFILE:-apptainer}" \
    ${CLI_OPTS:-""}

echo ""
echo "=== Creating FASTQ symlinks into data/raw-data/fastq/ ==="
# Downloaded filenames include an experiment accession prefix (e.g. SRX25813195_SRR30355320_1.fastq.gz).
# Only symlinks for files that were actually downloaded are created.
create_fastq_symlinks() {
    local samplesheet="$1"
    local fastq_base="../../data/raw-data/fastq"
    local source_fastq_dir
    # Get absolute path for symlink
    source_fastq_dir="$(realpath "../../data/source/fetchngs/fastq")"
    # Read csv entry
    while IFS=',' read -r _sample _run_id fastq_1 fastq_2; do
        for fastq_path in "${fastq_1}" "${fastq_2}"; do
            [ -z "${fastq_path}" ] && continue  # -z: true if string is empty | skips fastq if no value
            local csv_basename
            csv_basename="$(basename "${fastq_path}")"
            # Match against the downloaded file which may carry an experiment accession prefix
            local actual_file
            actual_file=$(find "${source_fastq_dir}" -maxdepth 1 -name "*${csv_basename}")
            [ -z "${actual_file}" ] && continue  # -z: true if string is empty | Don't make link if file isn't downloaded
            local dest="${fastq_base}/${fastq_path}"
            [ -e "${dest}" ] && continue  # -e: true if dest exists and is valid (false for broken symlinks, so those get recreated)
            mkdir -p "$(dirname "${dest}")"
            ln -sf "${actual_file}" "${dest}"
        done
    done < <(tail -n +2 "${samplesheet}")
}
create_fastq_symlinks "../../data/raw-data/fastq/low_pass_samples.csv"
create_fastq_symlinks "../../data/raw-data/fastq/high_pass_samples.csv"
echo "=== Symlinks complete ==="

echo ""
echo "=== Cleaning up obsolete work directories ==="
nextflow clean -f -before last && find work -type d -empty -delete
echo "=== Workflow complete ==="
