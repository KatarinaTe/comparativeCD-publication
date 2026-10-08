#! /usr/bin/env bash
#SBATCH -A uppmax2025-2-420
#SBATCH --nodes=1
#SBATCH --exclusive
# --mem must cover the largest resourceLimits.memory declared among the nested pipelines this
# job runs serially in-process (fetchngs' own nextflow.config apptainer profile: 750.GB) -- 16G
# here (fixed 2026-09-24 after a real 13h OOM kill on batch_05/low, sacct: ReqMem=16G,
# OUT_OF_ME+) silently capped the job far below what Nextflow's local executor believed it could
# schedule concurrently, even though --exclusive already grants the whole node (sinfo: 772000 MB
# on Pelle's p36-p49). 750G leaves ~22GB headroom below the node's real total, matching what the
# apptainer profile itself already assumes is safe.
#SBATCH --mem=750G
#SBATCH -t 7-00:00:00
#SBATCH -J fetch-raw-data-batched
#SBATCH --chdir=.

# Batched replacement for run_nextflow.sh's full-cohort fetch. run_nextflow.sh (and its
# params.yml) stay as-is and are still useful for small/manual single-batch fetches — this
# script is what actually gets `sbatch`-ed for a full 3,043-run cohort fetch+map.
#
# Why this exists: fetching all 3,043 SRA runs (11 HighPass + 3,032 LowPass) in one
# nf-core/fetchngs invocation, with 01_mapping's original (also-symlink-publishing)
# combined pipeline behind it, grows disk usage monotonically for the whole run with no
# safe way to reclaim space until everything finished — which never happened before a
# 15,000 GiB HPC project quota was exhausted. See archive/conversion_notes/01_mapping.md
# for the full writeup.
#
# This script instead, per track (HighPass: one batch of all 11 samples; LowPass: batches
# of BATCH_SIZE samples, see below) and per batch:
#   1. Slices that batch's rows out of {high,low}_pass_samples.csv and builds a matching
#      SRA accession sub-list (batch_state/<track>/batch_NN/).
#   2. Runs nf-core/fetchngs for just that batch's accessions.
#   3. Builds that batch's 1_fetch_map samplesheet from fetchngs' batch output (same
#      basename-matching logic as run_nextflow.sh's own create_fastq_symlinks()).
#   4. Runs code/01_mapping/1_fetch_map for that batch.
#   5. Materializes (dereference-copies) the batch's BAM/BAI/gVCF/gVCF_TBI/depth-stats
#      files into the permanent archive data/raw-data/mapped/<track>/ (one flat directory
#      per track, all samples' files side by side — every filename is already
#      sample-prefixed, so nothing collides — matching what 03_qc/4_qc_signals_IGV's
#      bam_dir params expect: a single directory to glob/lookup by sample_id, not a
#      per-sample subdirectory), and appends a row to the running
#      data/raw-data/mapped/<track>_pass_cohort_samplesheet.csv.
#   6. Deletes the batch's fetchngs and 1_fetch_map output/work directories — only after
#      step 5 is confirmed to have actually produced every expected file.
# This script's job ends once every batch (both tracks) is fetched+mapped+materialized — it does
# NOT also run code/01_mapping/2_cohort_genotyping (GENOTYPE_CALLING joint-calls all HighPass
# gVCFs together; GLIMPSE_LIKELIHOODS needs genotype likelihoods computed jointly across ALL
# LowPass BAMs at once — running it per-batch would silently produce different, less accurate
# results, not just run faster). That's a deliberately separate, manually-triggered step
# (`pixi run 01b-cohort-genotyping`), not auto-chained here as of 2026-10-05 -- it previously was,
# but: (1) given how often a batch needs a retry in practice, the auto-chain's "only run once
# literally everything succeeds in one job" condition rarely fires anyway; (2) cohort genotyping's
# own resource/walltime profile has no reason to match this job's (`--exclusive --mem=750G -t
# 7-00:00:00`, sized for fetchngs specifically); (3) a separate job is far easier to monitor and
# resubmit on its own terms than a step buried at the tail of a week-long fetch job. See
# `01b-cohort-genotyping` in pixi.toml.
#
# Resumable: each batch is marked done (batch_state/<track>/batch_NN.done) only after its outputs
# are materialized. A batch that exhausts its retries is skipped, not fatal -- the script moves on
# to the next one and reports every batch still needing a retry at the end (see FAILED_BATCHES
# below), rather than one bad batch stopping forward progress on every other batch in the same
# job's wall-clock budget. Re-submitting this script after a kill/timeout/crash (or just because
# some batches failed) skips every batch already marked done and continues from where it left off
# — it does not redo completed work. batch_state/ itself (small CSVs + marker files) is never
# deleted by this script.

set -euo pipefail

trap 'echo "=!= ERROR: run_batched.sh failed (exit code: $?) ==="' ERR

# Anchor to the submission directory — see run_nextflow.sh's own comment on why "$0"/
# BASH_SOURCE can't be used under sbatch.
cd "${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"

# ── Configuration ────────────────────────────────────────────────────────────

BATCH_SIZE="${BATCH_SIZE:-150}"          # LowPass samples per batch
HIGH_PASS_BATCH_SIZE=100000              # effectively unbounded: HighPass is always one batch (11 samples)
FETCHNGS_VERSION=1.12.0
# Profile for the inline `nextflow run nf-core/fetchngs` call below — this dir's own
# nextflow.config only defines `apptainer`/`test_one_sample`, same default as run_nextflow.sh.
FETCHNGS_PROFILE="${FETCHNGS_PROFILE:-apptainer}"
# Profile forwarded to both code/01_mapping/*/run_nextflow.sh calls — their own
# nextflow.configs default to `uppmax` (nf-core/configs' Rackham profile) when unset;
# kept as a separate default here so setting NXF_PROFILE for this script doesn't also
# have to know fetchngs' unrelated profile name.
NXF_PROFILE="${NXF_PROFILE:-uppmax}"

# Added 2026-10-02 after two separate, transient NCBI SRA resolver errors ("Cannot nest groups
# near index 72" -- a malformed-query response from `prefetch`'s internal SDL lookup, confirmed
# via ENA's own API and NCBI eutils to NOT reflect anything actually wrong with the accession
# itself) each killed the entire multi-day batched job, requiring a manual resubmit. See
# run_fetchngs_with_retry() below.
FETCHNGS_MAX_ATTEMPTS="${FETCHNGS_MAX_ATTEMPTS:-3}"            # 1 initial try + up to 2 retries
FETCHNGS_RETRY_DELAY_SECONDS="${FETCHNGS_RETRY_DELAY_SECONDS:-300}"

SRA_ACCESSION_LIST="./sra_accession_list.csv"
FASTQ_DIR="../../data/raw-data/fastq"
MAPPED_DIR="../../data/raw-data/mapped"
FETCHNGS_BATCH_ROOT="../../data/source/fetchngs_batches"
FETCH_MAP_BATCH_ROOT="../../data/results/01_mapping/1_fetch_map/batches"
FETCH_MAP_DIR="../01_mapping/1_fetch_map"
STATE_DIR="./batch_state"

# Filename suffixes produced by 1_fetch_map's DOG10K_MERGE_MARKDUPS_BQSR_GVCF/RUN_STATS
# (see code/01_mapping/1_fetch_map/process_definitions.nf) — samples are prefixed with
# their sample_id, e.g. <sample>.sorted.merged.MarkDups.BQSR.bam.
BAM_SUFFIX=".sorted.merged.MarkDups.BQSR.bam"
GVCF_SUFFIX=".sorted.merged.MarkDups.BQSR.g.vcf.gz"
DEPTH_SUFFIX=".sorted.merged.MarkDups.BQSR.bam.knownsites.depth.txt"

# ── Idempotency helpers ──────────────────────────────────────────────────────

is_done() { # track batch_id
    [ -e "${STATE_DIR}/$1/$2.done" ]
}

mark_done() { # track batch_id
    mkdir -p "${STATE_DIR}/$1"
    date -u +"%Y-%m-%dT%H:%M:%SZ" > "${STATE_DIR}/$1/$2.done"
}

abspath() { # path (must already exist)
    ( cd "$(dirname "$1")" && printf '%s/%s\n' "$(pwd)" "$(basename "$1")" )
}

# ── Batch splitting ──────────────────────────────────────────────────────────

# split_into_batches TRACK SAMPLESHEET BATCH_SIZE
# Slices SAMPLESHEET (header + `sample,run_id,fastq_1,fastq_2` rows, same format as
# data/raw-data/fastq/{high,low}_pass_samples.csv) into BATCH_SIZE-row batches under
# batch_state/TRACK/batch_NN/samplesheet.csv, plus a matching SRA accession sub-list
# (batch_state/TRACK/batch_NN/accessions.csv — one accession per line, no header, same
# format as sra_accession_list.csv) built from the samplesheet's `sample` column (col 1),
# filtered to only the accessions that are actually present in sra_accession_list.csv.
# A small number of LowPass rows carry no real SRA accession at all (see
# data/raw-data/fastq/labels_with_missing_sra_accession.txt — 12 rows, `sample` and
# `run_id` both set to a non-SRR Darwin's Ark ID); those are silently excluded from
# accessions.csv here, the same way create_fastq_symlinks() in run_nextflow.sh already
# silently skips any row whose FASTQ was never downloaded. Already-populated batch dirs
# are left untouched (idempotent — safe to call again on resume).
split_into_batches() {
    local track=$1 samplesheet=$2 batch_size=$3
    local track_state="${STATE_DIR}/${track}"
    mkdir -p "${track_state}"

    local header
    header=$(head -n1 "${samplesheet}")

    # NR>1: skip the header row | batch number = 1-indexed group of `batch_size` data rows
    awk -v n="${batch_size}" -v outdir="${track_state}" '
        NR>1 {
            batch = int((NR-2)/n) + 1
            bnum = sprintf("%02d", batch)
            print > (outdir "/.rows_" bnum)
        }
    ' "${samplesheet}"

    for rows_file in "${track_state}"/.rows_*; do
        [ -e "${rows_file}" ] || continue  # -e: true if glob matched a real file | no batches (empty samplesheet)
        local bnum bdir
        bnum=$(basename "${rows_file}" | sed 's/^\.rows_//')
        bdir="${track_state}/batch_${bnum}"
        mkdir -p "${bdir}"
        if [ ! -e "${bdir}/samplesheet.csv" ]; then
            { echo "${header}"; cat "${rows_file}"; } > "${bdir}/samplesheet.csv"
        fi
        if [ ! -e "${bdir}/accessions.csv" ]; then
            # -Fxf: fixed-string, whole-line match against patterns from sra_accession_list.csv
            # `|| true`: grep exits 1 (not an error here) when a batch has zero valid accessions
            cut -d',' -f1 "${rows_file}" | grep -Fxf "${SRA_ACCESSION_LIST}" - > "${bdir}/accessions.csv" || true
        fi
        rm -f "${rows_file}"
    done
}

# ── Per-batch samplesheet construction ───────────────────────────────────────

# build_fetch_map_samplesheet BATCH_SAMPLESHEET FETCHNGS_FASTQ_DIR OUT_CSV
# Mirrors create_fastq_symlinks() in run_nextflow.sh: fetchngs' downloaded FASTQ
# filenames carry an experiment-accession prefix (e.g.
# SRX25813195_SRR30355320_1.fastq.gz), matched here by basename suffix against each
# row's fastq_1/fastq_2 basename. FETCHNGS_FASTQ_DIR must be an absolute path so the
# resulting samplesheet resolves regardless of nf-schema's working directory. Rows
# whose FASTQ was never downloaded (not fetched by SRA, or one of the 12 samples with
# no real accession — see split_into_batches()) are skipped, not an error.
build_fetch_map_samplesheet() {
    local batch_samplesheet=$1 fetchngs_fastq_dir=$2 out_csv=$3
    echo "sample,run_id,fastq_1,fastq_2" > "${out_csv}"

    while IFS=',' read -r sample run_id fastq_1 fastq_2; do
        local f1_base f2_base f1_actual f2_actual
        f1_base="$(basename "${fastq_1}")"
        f2_base="$(basename "${fastq_2}")"
        f1_actual="$(find "${fetchngs_fastq_dir}" -maxdepth 1 -name "*${f1_base}" | head -n1)"
        f2_actual="$(find "${fetchngs_fastq_dir}" -maxdepth 1 -name "*${f2_base}" | head -n1)"
        [ -z "${f1_actual}" ] && continue  # -z: true if string is empty | R1 not downloaded, skip sample
        [ -z "${f2_actual}" ] && continue  # -z: true if string is empty | R2 not downloaded, skip sample
        # nf-schema's samplesheet CSV-to-JSON conversion infers a purely-numeric cell as a
        # JSON integer, which then fails schema_input.json's `run_id: {"type": "string"}`
        # check — common for LowPass rows, whose run_id is often a raw numeric Darwin's Ark
        # lab ID rather than a letter-prefixed SRR accession. Prefix to guarantee it never
        # looks like a number; safe because run_id is only ever used as a Nextflow
        # tag/read-group ID (code/01_mapping/1_fetch_map/process_definitions.nf), never
        # matched against external data.
        [[ "${run_id}" =~ ^[0-9]+$ ]] && run_id="run_${run_id}"
        echo "${sample},${run_id},${f1_actual},${f2_actual}" >> "${out_csv}"
    done < <(tail -n +2 "${batch_samplesheet}")
}

# ── Materialization ──────────────────────────────────────────────────────────

# materialize_batch TRACK FETCH_MAP_OUTDIR EXPECTED_COUNT
# Dereference-copies (cp -aL — verified on a placeholder symlink to actually copy real
# bytes, not another symlink, before relying on it here) this batch's published
# BAM/BAI/gVCF/gVCF_TBI/depth files out of FETCH_MAP_OUTDIR into the permanent archive
# data/raw-data/mapped/TRACK/ — one flat directory per track, every sample's files side
# by side (filenames are already sample-prefixed, so nothing collides across samples or
# batches) — matching what 03_qc/4_qc_signals_IGV's bam_dir params expect: a single
# directory to glob/look up by sample_id, not a per-sample subdirectory. Double-checks
# the copies are real files (not still symlinks), and appends one row per sample to
# data/raw-data/mapped/TRACK_pass_cohort_samplesheet.csv (created with its header if
# this is the first batch for TRACK). Returns non-zero — leaving the batch unmarked and
# its intermediates undeleted — if any expected file is missing or if fewer than
# EXPECTED_COUNT samples were found.
materialize_batch() {
    local track=$1 fetch_map_outdir=$2 expected_count=$3
    local mapped_track_dir="${MAPPED_DIR}/${track}"
    local cohort_csv="${MAPPED_DIR}/${track}_pass_cohort_samplesheet.csv"
    mkdir -p "${mapped_track_dir}"
    [ -e "${cohort_csv}" ] || echo "sample,bam,bai,gvcf,gvcf_tbi,depth_stats" > "${cohort_csv}"

    local produced=0
    while IFS= read -r bam_path; do
        local sample bai gvcf gvcf_tbi depth
        sample="$(basename "${bam_path}" "${BAM_SUFFIX}")"
        bai="${bam_path}.bai"
        gvcf="$(find "${fetch_map_outdir}/gvcf" -maxdepth 1 -name "${sample}${GVCF_SUFFIX}")"
        gvcf_tbi="${gvcf:+${gvcf}.tbi}"
        depth="$(find "${fetch_map_outdir}/depth" -maxdepth 1 -name "${sample}${DEPTH_SUFFIX}")"

        for f in "${bam_path}" "${bai}" "${gvcf}" "${gvcf_tbi}" "${depth}"; do
            if [ -z "${f}" ] || [ ! -e "${f}" ]; then
                echo "=!= ERROR: [${track}] expected output missing for sample ${sample} (looked for: '${f}') ===" >&2
                return 1
            fi
        done

        cp -aL "${bam_path}" "${bai}" "${gvcf}" "${gvcf_tbi}" "${depth}" "${mapped_track_dir}/"

        for f in "${bam_path}" "${bai}" "${gvcf}" "${gvcf_tbi}" "${depth}"; do
            local dest="${mapped_track_dir}/$(basename "${f}")"
            if [ -L "${dest}" ]; then  # -L: true if dest is a symlink (materialization failed to dereference)
                echo "=!= ERROR: [${track}] materialized file is still a symlink, not dereferenced: ${dest} ===" >&2
                return 1
            fi
        done

        # Idempotent on retry (e.g. after a crash mid-batch): a sample already recorded
        # from an earlier partial attempt at this same batch is not re-appended, even
        # though its files are harmlessly re-copied above.
        if ! grep -q "^${sample}," "${cohort_csv}"; then
            echo "${sample},$(abspath "${mapped_track_dir}/$(basename "${bam_path}")"),$(abspath "${mapped_track_dir}/$(basename "${bai}")"),$(abspath "${mapped_track_dir}/$(basename "${gvcf}")"),$(abspath "${mapped_track_dir}/$(basename "${gvcf_tbi}")"),$(abspath "${mapped_track_dir}/$(basename "${depth}")")" >> "${cohort_csv}"
        fi
        produced=$((produced+1))
    done < <(find "${fetch_map_outdir}/bam" -maxdepth 1 -name "*${BAM_SUFFIX}")

    if [ "${produced}" -lt "${expected_count}" ]; then
        echo "=!= ERROR: [${track}] only ${produced}/${expected_count} samples materialized ===" >&2
        return 1
    fi
}

# ── fetchngs retry wrapper ───────────────────────────────────────────────────

# run_fetchngs_with_retry ACCESSIONS OUTDIR WORKDIR
# Runs nf-core/fetchngs for one batch. Two layers of retry, added/revised 2026-10-05 after
# SRATOOLS_PREFETCH kept intermittently failing against NCBI's SRA Data Locator ("Cannot nest
# groups near index 72" -- a backend query-parser error, a different accession each time):
#
#   1. Task-level: fetchngs_overrides.config (-c'd in below) sets errorStrategy='retry',
#      maxRetries=5 on just the SRATOOLS_PREFETCH process. Handles the common case -- one
#      accession's fetch failing transiently -- by retrying only that task, leaving its ~149
#      unaffected siblings running. This is the main fix; the whole-invocation retry below is now
#      a backstop for whatever isn't absorbed at the task level.
#   2. Whole-invocation: up to FETCHNGS_MAX_ATTEMPTS retries of this entire function, with a fixed
#      FETCHNGS_RETRY_DELAY_SECONDS pause between them, same as before -- now with bare `-resume`
#      (see below) so a retry doesn't re-fetch accessions that already succeeded.
#
# Bare `-resume`, no explicit `-name`/target -- two real attempts at something more "deterministic"
# both turned out to be unnecessary and actively broken, worth recording so it isn't retried a
# third time:
#   - Originally omitted entirely, reasoned as -- `-resume` without an explicit target resumes the
#     *last recorded session in this directory* regardless of which batch it belonged to, since
#     this function launches from the same directory (analyses/00_fetch-raw-data) for every batch
#     in sequence -- too easy to accidentally resume a different batch's run.
#   - 2026-10-05 first fix: added a deterministic `-name "fetchngs_${track}_${batch_id}"` +
#     `-resume` to remove that ambiguity. Broke the job for real: job 7089448 died after 3 retries
#     all hitting `AbortOperationException: Run name 'fetchngs_low_batch_03' has been already
#     used` -- `-name` always means "create a new history entry with this name," and rejects it
#     outright if that name already exists, regardless of what `-resume` targets.
#   - 2026-10-05 same-day second attempt: tried `-resume` given that SAME name as its target
#     (instead of bare) -- still hit the identical collision error. `-resume <name>` does not make
#     an already-used `-name <name>` acceptable; they're independent concerns.
#   - Tested directly (local trivial pipeline, isolated -work-dir per "batch", genuinely differing
#     task inputs) before settling on this: bare `-resume`, issued right after a *different*
#     batch's run, correctly does NOT reuse that other batch's cached result (confirmed
#     `cached=0`) -- Nextflow's actual cache match is by task hash (inputs + script), not by which
#     session `-resume` happens to land on. The failure mode the original comment worried about
#     only reproduces when two tasks have *identical* inputs, which doesn't happen here: each
#     batch's `accessions.csv` is a disjoint slice of the full cohort, so no two batches' fetch
#     tasks ever hash the same. Re-running the *same* batch's command afterward correctly resumed
#     (`cached=1`) -- the actual behaviour this whole mechanism exists for.
#
# Only a failure that outlasts every whole-invocation attempt propagates via `set -e` and fails
# the job -- same behaviour as before this wrapper existed.
run_fetchngs_with_retry() {
    local accessions=$1 outdir=$2 workdir=$3
    local attempt=1
    while true; do
        echo "--- [nf-core/fetchngs] attempt ${attempt}/${FETCHNGS_MAX_ATTEMPTS} ---"
        # fetchngs 1.12.0's DSL2 predates Nextflow's strict syntax parser (now default) — same
        # pin as run_nextflow.sh's own fetchngs call. Scoped to just this command (not exported)
        # so it doesn't leak into 1_fetch_map's own run_nextflow.sh below.
        # set +e/-e bracket just this call so a failed attempt is handled by this loop instead
        # of immediately tripping the script's own ERR trap.
        set +e
        NXF_SYNTAX_PARSER=v1 nextflow run nf-core/fetchngs \
            -r "${FETCHNGS_VERSION}" \
            -params-file params.yml \
            -c "$(abspath fetchngs_overrides.config)" \
            --input "$(abspath "${accessions}")" \
            --outdir "${outdir}" \
            -work-dir "${workdir}" \
            -profile "${FETCHNGS_PROFILE}" \
            -resume \
            ${CLI_OPTS:-""}
        local status=$?
        set -e

        [ "${status}" -eq 0 ] && return 0
        if [ "${attempt}" -ge "${FETCHNGS_MAX_ATTEMPTS}" ]; then
            echo "=!= ERROR: nf-core/fetchngs failed after ${FETCHNGS_MAX_ATTEMPTS} attempts (exit ${status}) ===" >&2
            return "${status}"
        fi
        echo "=== nf-core/fetchngs attempt ${attempt} failed (exit ${status}) — retrying in ${FETCHNGS_RETRY_DELAY_SECONDS}s ==="
        sleep "${FETCHNGS_RETRY_DELAY_SECONDS}"
        attempt=$((attempt+1))
    done
}

# ── Per-batch driver ─────────────────────────────────────────────────────────

run_batch() { # track batch_id
    local track=$1 batch_id=$2
    local bdir="${STATE_DIR}/${track}/${batch_id}"
    local accessions="${bdir}/accessions.csv"
    local batch_samplesheet="${bdir}/samplesheet.csv"
    local fetch_map_samplesheet="${bdir}/fetch_map_samplesheet.csv"

    local fetchngs_batch_root="${FETCHNGS_BATCH_ROOT}/${track}/${batch_id}"
    local fetchngs_outdir="${fetchngs_batch_root}/output"
    local fetchngs_workdir="${fetchngs_batch_root}/work"
    local fm_batch_root="${FETCH_MAP_BATCH_ROOT}/${track}/${batch_id}"
    local fm_outdir="${fm_batch_root}/output"
    local fm_workdir="${fm_batch_root}/work"

    if [ ! -s "${accessions}" ]; then
        echo "=== [${track}/${batch_id}] No valid SRA accessions in this batch, nothing to fetch — marking done ==="
        mark_done "${track}" "${batch_id}"
        return
    fi

    echo "=== [${track}/${batch_id}] Fetching $(wc -l < "${accessions}") run(s) via nf-core/fetchngs ==="
    mkdir -p "${fetchngs_outdir}" "${fetchngs_workdir}"
    run_fetchngs_with_retry "${accessions}" "${fetchngs_outdir}" "${fetchngs_workdir}"

    echo "=== [${track}/${batch_id}] Building 1_fetch_map samplesheet ==="
    build_fetch_map_samplesheet "${batch_samplesheet}" "$(abspath "${fetchngs_outdir}/fastq")" "${fetch_map_samplesheet}"

    echo "=== [${track}/${batch_id}] Running 1_fetch_map ==="
    mkdir -p "${fm_outdir}" "${fm_workdir}"
    # Resolve abspaths *before* the subshell's cd below — abspath() resolves relative to
    # the current directory, and these paths are relative to this script's own cwd
    # (analyses/00_fetch-raw-data), not to FETCH_MAP_DIR.
    local abs_fetch_map_samplesheet abs_fm_workdir abs_fm_outdir
    abs_fetch_map_samplesheet="$(abspath "${fetch_map_samplesheet}")"
    abs_fm_workdir="$(abspath "${fm_workdir}")"
    abs_fm_outdir="$(abspath "${fm_outdir}")"
    (
        cd "${FETCH_MAP_DIR}"
        SAMPLESHEET="${abs_fetch_map_samplesheet}" \
        WORK_DIR="${abs_fm_workdir}" \
        OUTPUT_DIR="${abs_fm_outdir}" \
        NXF_PROFILE="${NXF_PROFILE}" \
        ./run_nextflow.sh
    )

    echo "=== [${track}/${batch_id}] Materializing outputs ==="
    materialize_batch "${track}" "${fm_outdir}" "$(wc -l < "${accessions}")"

    echo "=== [${track}/${batch_id}] Cleaning up batch intermediates ==="
    rm -rf "${fetchngs_batch_root}" "${fm_batch_root}"

    mark_done "${track}" "${batch_id}"
    echo "=== [${track}/${batch_id}] Done ==="
}

# ── Per-track driver ─────────────────────────────────────────────────────────

# Populated below when a batch exhausts its retries and is skipped rather than aborting the whole
# job -- checked by main() before running cohort genotyping, which needs the complete cohort, not
# just whatever got through in this job's wall-clock budget.
FAILED_BATCHES=()

process_track() { # track samplesheet batch_size
    local track=$1 samplesheet=$2 batch_size=$3
    echo "=== Splitting ${track} samplesheet (${samplesheet}) into batches of ${batch_size} ==="
    split_into_batches "${track}" "${samplesheet}" "${batch_size}"

    for bdir in "${STATE_DIR}/${track}"/batch_*; do
        [ -d "${bdir}" ] || continue
        local batch_id
        batch_id=$(basename "${bdir}")
        if is_done "${track}" "${batch_id}"; then
            echo "=== [${track}/${batch_id}] Already done, skipping ==="
            continue
        fi
        # run_batch in a subshell, not a bare call: `if ! run_batch ...; then` would disable
        # `set -e` for every command *inside* run_batch for the duration of this one call (a
        # classic bash gotcha -- errexit is suspended for an entire compound command used as an
        # if/while condition, not just the outermost one), which would let a genuine failure
        # partway through one batch's own steps (fetch succeeds, materialize fails, say) go
        # unnoticed instead of stopping that batch's processing at the right point. A subshell
        # keeps `set -e` fully active *inside* the attempt; only its overall exit status is
        # tested here, so a real failure still halts that one batch's own steps immediately, it
        # just doesn't take the rest of the job down with it.
        if ! ( run_batch "${track}" "${batch_id}" ); then
            echo "=!= WARNING: [${track}/${batch_id}] failed after exhausting retries -- skipping for" \
                 "now and continuing to the next batch. Re-submitting this job later retries only" \
                 "the batches that didn't complete; already-done ones (and any accessions already" \
                 "fetched within this one) are not redone. ===" >&2
            FAILED_BATCHES+=("${track}/${batch_id}")
            continue
        fi
    done
}

# ── Main ─────────────────────────────────────────────────────────────────────

main() {
    mkdir -p "${STATE_DIR}"

    process_track high "${FASTQ_DIR}/high_pass_samples.csv" "${HIGH_PASS_BATCH_SIZE}"
    process_track low  "${FASTQ_DIR}/low_pass_samples.csv"  "${BATCH_SIZE}"

    echo ""
    if [ "${#FAILED_BATCHES[@]}" -gt 0 ]; then
        echo "=== ${#FAILED_BATCHES[@]} batch(es) failed and were skipped:"
        printf '      %s\n' "${FAILED_BATCHES[@]}"
        echo "=== Re-submit this job to retry just the failed batch(es) -- already-done ones are"
        echo "=== skipped automatically, this does not restart from scratch."
        exit 1
    fi

    echo "=== All batches complete ==="
    echo "=== This script does not run cohort genotyping -- once every batch above is confirmed"
    echo "=== done, run it yourself: pixi run 01b-cohort-genotyping ==="
    echo "=== run_batched.sh complete ==="
}

# Guard so this file can be `source`d (e.g. by tests) to reuse the functions above
# without running the full batched fetch — mirrors the standard bash idiom for
# testable scripts. Under `sbatch`/direct execution, BASH_SOURCE[0] == $0 and main runs.
if [ "${BASH_SOURCE[0]}" = "${0}" ]; then
    main "$@"
fi
