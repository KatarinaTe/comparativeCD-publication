# Conversion notes — 01_mapping

Decisions and discrepancies noted while converting this stage's original scripts into the
Nextflow pipeline under `code/01_mapping/`. Moved out of `code/01_mapping/README.md` so that file
documents only the pipeline itself.

## Resolved Issues from Original Implementation

### RUN_STATS: duplication metrics crash

In the original `run-stats.py`, the duplication metrics file path is constructed with `.replace('.cram', '.sort.md.metricts.txt')` (note the typo `metricts`). The script is invoked with `.bam` files, so the replacement never matches and the expected metrics file path is never constructed. In the Nextflow `RUN_STATS` process the `dupMetFileName` construction and `summarize_stats()` call are commented out; all other metrics compute correctly.

### VQSR not applied

In the original `genomicsDB_joincall_SNPIDandHardFilter.sh`, `VariantRecalibrator` is called but `ApplyVQSR` is commented out. The pipeline falls through to hard filtering. This behaviour is preserved in `GENOTYPE_CALLING`.

### GATK version mismatch (original only)

The original `dog10k_mapping.sh` loads `module load GATK/4.1.4.1` but sets `gatk_path` to `gatk-4.2.0.0/`. Resolved in the Nextflow workflow by using a pinned container (GATK 4.4.0.0 via Seqera Wave).

### bcftools annotate output not compressed

In `genomicsDB_joincall_SNPIDandHardFilter.sh`, the `-Oz` compression flag is passed as a positional argument (`Oz`) rather than an option (`-Oz`). bcftools treats it as an input filename and writes uncompressed output. Fixed in `GENOTYPE_CALLING`.

### BCF → VCF conversion produces invalid output

In `04_LigatePerChromsome.sh`, `bcftools convert -O b` outputs BCF format (despite the misleading `.vcf` output filename), which is then passed through `bgzip`, producing a double-encoded file that bcftools cannot read as VCF. Fixed in `BCF_VCF` using `bcftools view -Oz`.

## Splitting into two pipelines (disk-quota crisis)

The original single `01_mapping.nf` (and `process_definitions.nf`) covered the entire stage —
FASTQ fetch/mapping through GLIMPSE imputation and joint SNP calling — as one Nextflow
invocation, with `publish_dir_mode: symlink` throughout. Combined with `analyses/00_fetch-raw-data`
also publishing `nf-core/fetchngs`'s FASTQ output as symlinks (never copying the downloaded
`.fastq.gz` out of `work/`), running the full 3,043-run cohort (11 HighPass + 3,032 LowPass, plus
12 LowPass samples with no real SRA accession — see `data/raw-data/fastq/labels_with_missing_sra_accession.txt`)
in one pass meant disk usage under `work/` grew monotonically for the run's entire duration, with
no safe point to reclaim space until everything finished. It never finished: a 15,000 GiB HPC
project quota was exhausted partway through. `nf-core/fetchngs`'s own SRA path independently makes
this worse — the `.sra` file, the uncompressed intermediate `.fastq`, and the final `.fastq.gz` can
all exist in `work/` simultaneously per run, on top of no concurrency cap (`executor.name = 'local'`).

The fix keeps `nf-core/fetchngs` as an external, unmodified pipeline (no vendoring of SRA-download
logic) and batches the whole fetch+map process: fetch a batch's FASTQs, map that batch, materialize
(dereference-copy) its outputs to a permanent location, delete the batch's intermediates, repeat.
This requires the mapping stage itself to be splittable, which is why `01_mapping.nf` was split into
two pipelines instead of just wrapping the existing one in a batching loop.

### Where the split was drawn

The dividing line is exactly where the original pipeline's data dependencies stop being
per-sample-independent:

- **`1_fetch_map/`** — `BWA_MEM2_INDEX`, `DOG10K_MAPPING`, `DOG10K_MERGE_MARKDUPS_BQSR_GVCF`,
  `RUN_STATS`. Every one of these operates on a single sample (or a single run, before the
  per-sample merge) in isolation — nothing here reads another sample's data. Safe to invoke any
  number of times over disjoint sample subsets, in any order, including re-invoking a batch that
  partially failed (each `run_nextflow.sh` sets `resume = true` and a batch-specific `-work-dir`,
  so a re-run resumes rather than redoing completed tasks).
- **`2_cohort_genotyping/`** — `BWA_MEM2_INDEX` (duplicate, see below), `GENOTYPE_CALLING`,
  `EXTRACT_DEPTH`, `GLIMPSE_DEFINE_CHUNKS`, `GLIMPSE_VARIABLE_SITES`, `GLIMPSE_LIKELIHOODS`,
  `GLIMPSE_IMPUTE_CHUNKS`, `INDEX_BCFS`, `GLIMPSE_LIGATE`, `BCF_VCF`, `GLIMPSE_QC_FILTER`,
  `VCF_TO_PLINK`, `MERGE_CHR_PLINK`. Everything here is cohort-wide by construction:
  - `GENOTYPE_CALLING` joint-calls all HighPass gVCFs together (GenomicsDB import + joint
    genotyping) — splitting the input gVCF set across batches would joint-call a different,
    smaller cohort each time, silently changing the called genotypes, not just deferring work.
  - `GLIMPSE_LIKELIHOODS` computes genotype likelihoods from the bundle of *all* LowPass BAMs at
    once, per chromosome. This is a hard correctness constraint, not a performance one: GLIMPSE's
    imputation accuracy depends on genotype likelihoods being computed jointly across the whole
    reference cohort — running this over a subset of samples would silently produce different,
    less accurate imputed genotypes for every sample in that subset, with no error or warning.

  Consequently `2_cohort_genotyping` must run exactly once, only after every `1_fetch_map` batch
  (both tracks) has completed and been materialized. Running it more than once, or before the
  cohort is complete, does not fail loudly — it just produces silently-wrong joint-calling/
  imputation results, which is the reason `analyses/00_fetch-raw-data/run_batched.sh` gates it
  behind an explicit `cohort_genotyping.done` marker.

  `BWA_MEM2_INDEX` is duplicated (byte-identical) into `2_cohort_genotyping/process_definitions.nf`
  purely so `2_cohort_genotyping.nf` can rebuild the same `ch_assembly_ref` tuple shape
  (`[fasta, [fai, dict, bwa-mem2 index files...]]`) that every downstream process signature still
  expects, unchanged from the original combined pipeline. No process in `2_cohort_genotyping`
  actually uses the bwa-mem2 index files themselves (GATK/bcftools/GLIMPSE only touch `.fai`/
  `.dict`) — this is flagged in the process's own header comment as redundant-but-harmless,
  kept for consistency rather than reworking every downstream process's input tuple shape.

Total processes: 16 in the original file, 16 unique across the split (4 unique to
`1_fetch_map` + 12 unique to `2_cohort_genotyping`, `BWA_MEM2_INDEX` shared/duplicated) —
confirmed nothing was lost or added in the split; every process's script body (`script:`/`stub:`
blocks) was diffed against the original file and is byte-identical.

### Symlink publishing and materialization

Both split pipelines keep `publish_dir_mode: symlink` (the default) exactly as the original
combined pipeline did — this is unchanged and still the right choice *within* a single pipeline
run, since it avoids doubling disk usage between `work/` and the published output during that
run. The crisis was never symlink-publishing itself; it was that nothing downstream of
`1_fetch_map` ever converted those symlinks into real, permanent files, so `work/` (and the
matching `fetchngs` `work/`) could never be deleted for a batch that had already finished.

`analyses/00_fetch-raw-data/run_batched.sh` is what closes that gap: after each batch's
`1_fetch_map` run completes, it dereference-copies (`cp -aL` — verified against a placeholder
symlink to confirm it copies real bytes, not another symlink, before relying on it for real BAM/
gVCF files) the batch's `bam/`, `gvcf/`, and `depth/` output into
`data/raw-data/mapped/{high,low}/` (one flat directory per track — every sample's files side by
side, since filenames are already sample-prefixed; fixed 2026-09-15 from an earlier per-sample
`<sample>/` subdirectory layout that didn't match what `03_qc`/`4_qc_signals_IGV` actually expect),
a permanent, non-symlinked archive. Only once those
copies are confirmed present (and one row per sample has been appended to that track's running
`data/raw-data/mapped/{high,low}_pass_cohort_samplesheet.csv`) does it delete that batch's
`fetchngs` and `1_fetch_map` output/work directories. `2_cohort_genotyping`'s two cohort
samplesheets (`assets/schemas/schema_cohort_input.json`) are how it reads those materialized,
already-merged per-sample files back in as its own input, instead of consuming an in-memory
channel from the same workflow run the way the original combined pipeline did.
