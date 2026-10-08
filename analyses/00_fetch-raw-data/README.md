# 00. Fetch Raw Data

Downloads FASTQ files from SRA and reference files from public repositories.

## Files

| File | Description |
|---|---|
| [`run_nextflow.sh`](run_nextflow.sh) | Runs nf-core/fetchngs to download FASTQ files, then creates symlinks into `data/raw-data/fastq/` — for small/manual fetches only, see "Disk space required" below |
| [`run_batched.sh`](run_batched.sh) | Batched fetch+map driver for a full-cohort run — see "Disk space required" below |
| [`fetchngs_overrides.config`](fetchngs_overrides.config) | Extra Nextflow config `-c`'d into every `run_batched.sh` fetchngs call — see "Reliability: NCBI fetch failures" below |
| [`fetch-reference-data.sh`](fetch-reference-data.sh) | Downloads reference files into `data/raw-data/references/`, and derives `data/raw-data/genes6_UU_Cfam_GSD_1.0_ROSY.txt` (GWAS clumping gene ranges) from the Dog10K gene models |
| [`fetch-scilifelab-deposit.sh`](fetch-scilifelab-deposit.sh) | Downloads our SciLifeLab deposit's pipeline inputs (doi.org/10.17044/scilifelab.33339309): Axiom genotypes into `data/raw-data/axiom/`, survey responses and breed/sex metadata into `data/02_data_processing/` — `pixi run 00c-fetch-scilifelab-deposit` |
| [`params.yml`](params.yml) | nf-core/fetchngs parameters (input accession list, output directory) |
| [`nextflow.config`](nextflow.config) | Executor profiles and resource limits |
| [`sra_accession_list.csv`](sra_accession_list.csv) | SRA accession numbers for all samples |

## Running

For a full cohort run (the normal case — see "Disk space required" below for why this is
different from what you might expect):

```bash
pixi run 00a-fetch-reference-data  # reference files only — small, run this first
pixi run 00b-mapping-batched    # fetches + maps the full 3,043-run cohort, in disk-bounded batches
```

`mapping-batched` (`run_batched.sh`) submits itself as a full-node SLURM batch job (14-day wall
time) and is fire-and-forget — the command returns as soon as the job is queued. Monitor progress
with `squeue -u $USER` or by tailing the `slurm-<jobid>.out` file written into this directory. It
fetches sequence data itself, per batch — there is no separate "fetch all the FASTQ first" step
for a full run. See "Disk space required" below for what it does and why.

`pixi run manual-fetch-sequence-data` (`run_nextflow.sh`, fetches ALL sequence data in one
`nf-core/fetchngs` run, symlinking FASTQ into `data/raw-data/fastq/{high,low}/`) still exists, but
**only use it for small/manual fetches** — e.g. a handful of samples for local testing. Using it
for the full cohort is what originally exhausted a 15,000 GiB HPC project quota; see below.

## Disk space required

`sra_accession_list.csv` has 3,043 runs: 11 HighPass (`data/raw-data/fastq/high_pass_samples.csv`)
and 3,032 LowPass (`data/raw-data/fastq/low_pass_samples.csv`). Based on each run's own SRA-reported
size (`total_size`, via NCBI's `esummary`, a reasonable proxy for download size — not necessarily
identical to the final `.fastq.gz` size fetchngs produces):

| | Runs | Measured/sampled | Estimated total |
|---|---|---|---|
| HighPass | 11 | All 11 queried directly | **~680 GB** |
| LowPass | 3,032 | 15 runs sampled at random (mean ~0.96 GB/run) | **~2.9 TB** (extrapolated, not exhaustive) |
| **Total** | 3,043 | | **~3.5 TB** |

This is download size only. `nf-core/fetchngs` downloads the SRA-format file first, then converts
it to `.fastq.gz` — during that conversion, both the source file and its FASTQ output can exist on
disk at once, so **peak working-directory usage can exceed the ~3.5 TB total above** unless
intermediates are cleaned up as each run completes (fetchngs does this automatically per-run by
default; don't disable that cleanup to "save time" or peak usage will approach 2x).

**To avoid exhausting a shared project quota, use [`run_batched.sh`](run_batched.sh) instead of
`run_nextflow.sh` for a full-cohort fetch:**

```bash
sbatch run_batched.sh
```

It fetches + maps HighPass (one batch, all 11 samples) and LowPass (batches of `BATCH_SIZE`
samples, default 150 — override by setting the `BATCH_SIZE` environment variable before
submitting) sequentially, materializing and cleaning up each batch's intermediates before
starting the next, so disk usage stays bounded to roughly one batch's footprint at a time
instead of the full ~3.5 TB. A batch that exhausts its retries (see "Reliability" below) is
skipped, not fatal — the script keeps going on the remaining batches and reports every batch
still needing a retry at the end, rather than one bad batch stopping all forward progress for the
rest of the job's wall-clock budget. It is resumable — re-submitting (whether after a
kill/timeout, or just to retry batches that failed) picks up from the last completed batch rather
than redoing finished work.

**It does not run cohort genotyping.** `code/01_mapping/2_cohort_genotyping` is a deliberately
separate, manually-triggered step (`pixi run 01b-cohort-genotyping`) — run it yourself once every
batch above is confirmed done. It was auto-chained here until 2026-10-05; removed because (1) given
how often a batch needs a retry in practice, the old "only run once literally everything succeeds
in one job" condition rarely fired anyway, (2) cohort genotyping's own resource/walltime needs have
no reason to match this job's (sized for fetchngs specifically), and (3) a separate job is far
easier to monitor and resubmit on its own terms than a step buried at the tail of a week-long fetch
job.

See [`code/01_mapping/README.md`](../../code/01_mapping/README.md) and
[`archive/conversion_notes/01_mapping.md`](../../archive/conversion_notes/01_mapping.md) for why
this exists and how it's built.

`run_nextflow.sh`/`params.yml` are unchanged and still useful for small/manual fetches (e.g. just
the reference-adjacent data) — just not for a full-cohort run.

- The 11 HighPass runs alone are ~680 GB (roughly 20% of the total from 0.4% of the runs) — this
  is why HighPass and LowPass are always fetched as separate batches.
- Request quota sized to the ~3.5 TB estimate above (with headroom for the peak-usage caveat), not
  to the number of runs — the two aren't proportional here. `run_batched.sh` bounds *working*
  disk usage, not the eventual size of the permanent `data/raw-data/mapped/` archive, which still
  ends up holding the full cohort's BAM/gVCF/depth files.

This estimate hasn't been refined beyond a 15-run sample for LowPass (NCBI's `esummary` API wasn't
reliably reachable for a larger batch when this was measured) — treat ~2.9 TB as an order-of-magnitude
figure, not an exact one.

## Reliability: NCBI fetch failures

`nf-core/fetchngs`'s `SRATOOLS_PREFETCH` process intermittently fails against NCBI's SRA Data
Locator service with `Cannot nest groups near index 72` — a backend query-parser error, hit on a
different accession each time, not specific to any one SRA run. Community reports of related
`prefetch` reliability problems under high concurrency (e.g.
[ncbi/sra-tools#560](https://github.com/ncbi/sra-tools/issues/560)) are consistent with this being
a load-related backend issue rather than anything wrong with our accessions or pipeline. Also
tried and ruled out: `download_method: ftp` (bypasses NCBI entirely via ENA's mirror) — the
`SRA_FASTQ_FTP` process's `wget` container has no `/etc/resolv.conf` at all, so it cannot resolve
any hostname, a deterministic 100% failure rather than a fix (see `params.yml`'s own comment for
the full download-method history).

Two layers of retry in `run_batched.sh`, as of 2026-10-05:

1. **Task-level** (`fetchngs_overrides.config`, `-c`'d into every fetchngs call): retries just the
   one failing `SRATOOLS_PREFETCH` task up to 5 times, leaving its ~149 unaffected sibling fetches
   in that batch running. Handles the common case directly, rather than needing the whole-batch
   retry below.
2. **Whole-invocation** (`run_fetchngs_with_retry` in `run_batched.sh`): up to
   `FETCHNGS_MAX_ATTEMPTS` (default 3) retries of the entire ~150-accession fetchngs invocation,
   now resumed (bare `-resume`, no explicit name/target — see that function's own header comment
   for why an explicit `-name` was tried twice and broke the job both times) rather than restarted
   from scratch — already-completed fetches are skipped. This is now a backstop for whatever the
   task-level retry doesn't absorb, not the primary mechanism.

If a batch still exhausts both layers and the job fails, re-submitting (`pixi run
00b-mapping-batched` or `sbatch run_batched.sh`) picks up from the last completed batch — per-batch
idempotency (`batch_state/<track>/<batch_id>.done`) is unaffected by any of this.
