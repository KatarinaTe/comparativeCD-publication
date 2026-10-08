# 01. Mapping

Maps raw sequencing reads into analysis-ready genotype data via two parallel tracks (HighPass/LowPass) that converge in [`code/02_data_processing/3_Merging_filtering`](../02_data_processing/3_Merging_filtering).

**Pipeline engine:** Nextflow (two separate pipelines, see below) **Containerized:** Yes (Seqera Wave) **Estimated runtime:** \[FILL IN\]

------------------------------------------------------------------------

## Two pipelines, not one

This stage is split into two independently-invoked Nextflow pipelines:

| Pipeline | Directory | Runs |
|------------------------|------------------------|------------------------|
| Fetch + map | [`1_fetch_map/`](1_fetch_map/) | Once per batch of samples, any number of times |
| Cohort genotyping | [`2_cohort_genotyping/`](2_cohort_genotyping/) | Exactly once, over the complete cohort |

`1_fetch_map` is per-sample-independent (FASTQ → mapped/merged/BQSR'd BAM + gVCF + depth stats) and safe to run repeatedly over disjoint batches of samples — HighPass and LowPass samples go through identical processes here, track is just which samplesheet a given invocation is pointed at. `2_cohort_genotyping` joint-calls all HighPass gVCFs together and computes LowPass genotype likelihoods jointly across every LowPass BAM at once (GLIMPSE imputation accuracy depends on this), so it must only ever run once, after every `1_fetch_map` batch has completed for both tracks.

[`analyses/00_fetch-raw-data/run_batched.sh`](../../analyses/00_fetch-raw-data/run_batched.sh) is the driver that fetches raw data in disk-bounded batches, runs `1_fetch_map` per batch, materializes each batch's output into a permanent archive, and runs `2_cohort_genotyping` once at the end. See [`archive/conversion_notes/01_mapping.md`](../../archive/conversion_notes/01_mapping.md) for why this split exists (a disk-quota crisis fetching + mapping the full cohort in one pass) and exactly how the process/channel split was drawn.

------------------------------------------------------------------------

## Pipeline Overview

The diagram below shows the full logical data flow across both pipelines — everything up to and including `RUN_STATS` runs inside `1_fetch_map` (once per batch); everything from `GENOTYPE_CALLING`/`EXTRACT_DEPTH` onward runs inside `2_cohort_genotyping` (once, over the accumulated cohort).

``` mermaid
%%{init: {'theme': 'neutral'}}%%
flowchart TD
    subgraph inputs["Inputs"]
        HP_IN(["11 HighPass FASTQs"])
        LP_IN(["3044 LowPass FASTQs"])
    end

    REF(["Reference Genome"])
    BQSR(["Known Variants"])
    GLIMPSE_REF(["Dog10K Phased Panel"])

    subgraph hp["Track 1 — HighPass (11 samples)"]
        direction TB
        HP_MAP["Map · bwa-mem2"]
        HP_MERGE["Merge · MarkDups · BQSR · gVCF"]
        HP_STATS["Quality Stats"]
        HP_JC["Joint Calling · GenomicsDB"]
        HP_FILTER["SNP Hard Filter"]
        HP_MAP --> HP_MERGE --> HP_STATS
        HP_MERGE --> HP_JC --> HP_FILTER
    end

    subgraph lp["Track 2 — LowPass (3044 samples)"]
        direction TB
        LP_MAP["Map · bwa-mem2"]
        LP_MERGE["Merge · MarkDups · BQSR · gVCF"]
        LP_STATS["Quality Stats"]
        LP_DEF["Define Imputation Chunks"]
        LP_SITES["Extract Variable Sites"]
        LP_GL["Genotype Likelihoods"]
        LP_IMPUTE["GLIMPSE Phase & Impute"]
        LP_INDEX["Index BCF"]
        LP_LIGATE["Ligate Chunks"]
        LP_BCF2VCF["BCF → VCF"]
        LP_QC["QC Filter · INFO ≥ 0.8"]
        LP_PLINK["VCF → PLINK"]
        LP_MERGEALL["Merge Chromosomes"]
        LP_MAP --> LP_MERGE --> LP_STATS
        LP_MERGE --> LP_GL
        LP_SITES --> LP_GL
        LP_GL --> LP_IMPUTE
        LP_DEF --> LP_IMPUTE
        LP_IMPUTE --> LP_INDEX --> LP_LIGATE
        LP_LIGATE --> LP_BCF2VCF --> LP_QC
        LP_QC --> LP_PLINK --> LP_MERGEALL
    end

    HP_IN --> HP_MAP
    LP_IN --> LP_MAP

    REF -.-> HP_MAP
    REF -.-> LP_MAP
    REF -.-> LP_GL
    BQSR -.-> HP_MERGE
    BQSR -.-> LP_MERGE
    GLIMPSE_REF -.-> LP_DEF
    GLIMPSE_REF -.-> LP_SITES
    GLIMPSE_REF -.-> LP_IMPUTE

    HP_OUT[("SNP.HF.ann.id.vcf.gz")]
    LP_OUT[("DA_IMP_GENCOVE_ALLCHR .bed/.bim/.fam")]

    HP_FILTER --> HP_OUT
    LP_MERGEALL --> LP_OUT

    NEXT(["02_data_processing"])
    HP_OUT --> NEXT
    LP_OUT --> NEXT
```

> **HighPass and LowPass tracks are independent and run in parallel.** Both tracks share the reference genome, BQSR sites, depth sites, and chunk definitions. The LowPass track additionally requires the Dog10K phased imputation panel.

------------------------------------------------------------------------

## Input Data

### Samplesheets

`1_fetch_map` and `2_cohort_genotyping` take different samplesheets — one per run (pre-mapping) and one per sample (post-mapping) respectively.

**`1_fetch_map`** — schema: [`1_fetch_map/assets/schemas/schema_input.json`](1_fetch_map/assets/schemas/schema_input.json)

| Column    | Description                                                     |
|----------------------------|--------------------------------------------|
| `sample`  | Sample identifier (SRA run accession)                           |
| `run_id`  | Run identifier (multiple rows per sample for multi-run samples) |
| `fastq_1` | Path to R1 FASTQ (relative to samplesheet)                      |
| `fastq_2` | Path to R2 FASTQ (relative to samplesheet)                      |

| Samplesheet | Samples |
|------------------------------------------|------------------------------|
| `data/raw-data/fastq/high_pass_samples.csv` | 11 HighPass (≥6× coverage) — run as one batch |
| `data/raw-data/fastq/low_pass_samples.csv` | 3044 LowPass (0.1–2× coverage) — run in batches, see `run_batched.sh`'s `BATCH_SIZE` |

> 12 LowPass rows carry no real SRA accession and are excluded — not available in SRA. See [`data/raw-data/fastq/labels_with_missing_sra_accession.txt`](../../data/raw-data/fastq/labels_with_missing_sra_accession.txt). `run_batched.sh` filters these out of each batch's accession list automatically; they never reach `1_fetch_map` and are excluded downstream the same way they always have been.

**`2_cohort_genotyping`** — schema: [`2_cohort_genotyping/assets/schemas/schema_cohort_input.json`](2_cohort_genotyping/assets/schemas/schema_cohort_input.json)

| Column | Description |
|----------------------------|--------------------------------------------|
| `sample` | Sample identifier |
| `bam` / `bai` | Merged, MarkDups/BQSR'd BAM + index (`1_fetch_map` output, materialized) |
| `gvcf` / `gvcf_tbi` | Per-sample gVCF + index (`1_fetch_map` output, materialized) |
| `depth_stats` | Per-sample depth-stats file (`1_fetch_map` output, materialized) |

| Samplesheet | Samples |
|------------------------------------------|------------------------------|
| `data/raw-data/mapped/high_pass_cohort_samplesheet.csv` | All 11 HighPass, one row each |
| `data/raw-data/mapped/low_pass_cohort_samplesheet.csv` | All LowPass, one row each |

Both cohort samplesheets are built incrementally by `run_batched.sh` — one row is appended per sample as each batch's `1_fetch_map` output is materialized (see Output Data below). They are not hand-built and do not exist until at least one batch has completed.

### Reference Files

| Parameter | File | Source |
|-----------------|-----------------------------------|---------------------|
| `assembly_ref` | `UU_Cfam_GSD_1.0_ROSY.fa` + `.fai`, `.dict`, bwa-mem2 index | https://github.com/jmkidd/dogmap |
| `known_variants` | `UU_Cfam_GSD_1.0.BQSR.DB.bed.gz` + `.tbi` | https://github.com/jmkidd/dogmap |
| `depth_sites` | `SRZ189891_722g.simp.header.CanineHD.names.GSD_1.0.filter.vcf.gz` + `.tbi` | https://github.com/jmkidd/dogmap |
| `chunks_dir` | `chunks_10/` (chunk_1.bed … chunk_10.bed) | https://github.com/jmkidd/dogmap |
| `glimpse_ref_panel` | `AutoAndXPAR.Dog10K.Phased.bcf` + `.csi` | https://doi.org/10.5281/zenodo.8084059 |

------------------------------------------------------------------------

## Output Data

### `1_fetch_map` (per batch, before materialization)

Published under whatever `OUTPUT_DIR` a given batch invocation used (symlinks into that batch's `work/`, matching the original combined pipeline's publish mode):

| Path | Description |
|-----------------------------|-------------------------------------------|
| `bam/<sample>.sorted.merged.MarkDups.BQSR.bam(.bai)` | Per-sample merged, MarkDups/BQSR'd BAM |
| `gvcf/<sample>.sorted.merged.MarkDups.BQSR.g.vcf.gz(.tbi)` | Per-sample gVCF |
| `depth/<sample>....knownsites.depth.txt` | Per-sample autosome + chrX coverage |
| `pipeline_info/versions.yml` | Per-process tool versions |

`run_batched.sh` dereference-copies these into the permanent, disk-safe archive `data/raw-data/mapped/{high,low}/<sample>/` after each batch, and appends one row per sample to `data/raw-data/mapped/{high,low}_pass_cohort_samplesheet.csv` — see Input Data above.

### `2_cohort_genotyping` (once, final)

| Path | Description |
|-------------------------|-----------------------------------------------|
| `highpass/snp_vcf/SNP.HF.ann.id.vcf.gz` | Jointly called, hard-filtered, annotated SNP VCF (11 samples) |
| `lowpass/plink/DA_IMP_GENCOVE_ALLCHR.{bed,bim,fam}` | Merged PLINK dataset (3044 samples, chr1–38, post-imputation) |
| `lowpass/vcf/DA.${chr}.merged.qc_modi.vcf.gz` | Per-chromosome QC-filtered, imputed low-pass cohort VCF (pre-PLINK-conversion; consumed by `03_qc` concordance studies) |
| `lowpass/depth/combined_knownsites.depth.txt` | Cohort-level per-sample autosome + chrX coverage |
| `pipeline_info/versions.yml` | Per-process tool versions |

Per-sample BAMs are not re-published here — they already live permanently under `data/raw-data/mapped/{high,low}/<sample>/` (see above), upstream of this pipeline.

------------------------------------------------------------------------

## Process Map

| \# | Process | Pipeline | Track | Description |
|-------|--------------|--------------|--------------|-------------------|
| 1 | `BWA_MEM2_INDEX` | `1_fetch_map` | Both | Index reference genome |
| 2 | `DOG10K_MAPPING` | `1_fetch_map` | Both | bwa-mem2 → sorted BAM |
| 3 | `DOG10K_MERGE_MARKDUPS_BQSR_GVCF` | `1_fetch_map` | Both | Merge runs → MarkDuplicates → BQSR → gVCF |
| 4 | `RUN_STATS` | `1_fetch_map` | Both | Quality metrics (flagstat, insert size, depth) |
| — | `BWA_MEM2_INDEX` (duplicate) | `2_cohort_genotyping` | — | Re-run for tuple-shape consistency only — not otherwise used by anything downstream, see the process's own header comment |
| 5 | `EXTRACT_DEPTH` | `2_cohort_genotyping` | LowPass | Cohort-level depth aggregation |
| 6 | `GENOTYPE_CALLING` | `2_cohort_genotyping` | HighPass | Joint calling + SNP hard filtering |
| 7 | `GLIMPSE_DEFINE_CHUNKS` | `2_cohort_genotyping` | LowPass | Define imputation chunks from reference panel |
| 8 | `GLIMPSE_VARIABLE_SITES` | `2_cohort_genotyping` | LowPass | Extract biallelic SNP sites |
| 9 | `GLIMPSE_LIKELIHOODS` | `2_cohort_genotyping` | LowPass | Compute genotype likelihoods (jointly, across all LowPass BAMs) |
| 10 | `GLIMPSE_IMPUTE_CHUNKS` | `2_cohort_genotyping` | LowPass | Phase and impute each chunk |
| 11 | `INDEX_BCFS` | `2_cohort_genotyping` | LowPass | Index imputed BCF files |
| 12 | `GLIMPSE_LIGATE` | `2_cohort_genotyping` | LowPass | Ligate chunks per chromosome |
| 13 | `BCF_VCF` | `2_cohort_genotyping` | LowPass | Convert BCF → bgzipped VCF |
| 14 | `GLIMPSE_QC_FILTER` | `2_cohort_genotyping` | LowPass | Filter on INFO score (≥ 0.8) |
| 15 | `VCF_TO_PLINK` | `2_cohort_genotyping` | LowPass | Convert VCF → PLINK binary |
| 16 | `MERGE_CHR_PLINK` | `2_cohort_genotyping` | LowPass | Merge all chromosomes |

16 unique processes total (`BWA_MEM2_INDEX` is duplicated across both pipeline directories, byte-identical) — matches the original combined `01_mapping.nf`'s process count exactly.

------------------------------------------------------------------------

## Known characteristics

- `RUN_STATS` does not compute duplication-rate metrics; all other depth/coverage metrics are computed as normal.
- Variant recalibration (VQSR) is not applied — `GENOTYPE_CALLING` uses hard-filtering only.
- GATK 4.4.0.0 is used throughout, via a pinned container.
- **No working stub test currently exists for either substage** (verified 2026-09-18: `pixi run 01a-fetch-map -stub` fails on a missing `SAMPLESHEET`; `pixi run 01b-cohort-genotyping -stub` fails schema validation on 4 missing input paths). Earlier fixtures this stage's stub test once relied on are gone, likely orphaned when the `run_batched.sh` disk-quota rework split this stage into `1_fetch_map`/`2_cohort_genotyping` and moved sample-sheet construction into that driver script. Not rebuilt, since real-data validation happens on the cluster instead — a fresh contributor relying on a stub test to sanity-check this stage without cluster access will hit this.
- **`DOG10K_MAPPING` can intermittently fail with `Shringking not supported yet` / `munmap_chunk(): invalid pointer`.** This is a confirmed, known upstream race condition in `bwa-mem2==2.1`'s multi-threaded memory-reallocation code ([bwa-mem2/bwa-mem2#88](https://github.com/bwa-mem2/bwa-mem2/issues/88)) — not a bug in this repo's code or your data, and not fixed in any released version since (pinned to match the manuscript's own stated version). It's genuinely probabilistic: the same input can succeed or fail on different runs, so `1_fetch_map/conf/base.config` auto-retries this specific process on exit status 1 (up to 3 attempts) — a real occurrence of this is logged in `REPRODUCIBILITY_AUDIT.md`'s 2026-09-23 update. If it still fails after 3 attempts on the same task, that's worth investigating further rather than just retrying again. Retrying alone wasn't reliably clearing the race on one sample that hit it on 3 separate retries in each of two independent full-cohort runs (2026-09-27, 2026-09-28), so `DOG10K_MAPPING`'s `cpus` was also dropped from 12 (the default `process_high` label) to 8 — narrowing the timing window while keeping most of the throughput — applied to every task in this process, not just the sample that triggered it.

------------------------------------------------------------------------

## Container

Five Seqera Wave containers, shared across both pipeline directories where a process appears in both (`BWA_MEM2_INDEX`): `dog10k_mapping` (bwa-mem2/samtools/picard/gatk/tabix — mapping, merge/BQSR/gVCF, quality stats), `gawk`, `highpass_snpfiltering`, `glimpse-bio`, `plink`. See each process's own `container` directive in `1_fetch_map/process_definitions.nf` / `2_cohort_genotyping/process_definitions.nf` for exact tags.