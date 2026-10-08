# 02. Data processing — 2. Axiom imputation

Lifts the Darwin's Ark Axiom SNP array genotypes from canFam3 to canFam4, then imputes each chromosome against the Dog10K phased reference panel, producing a merged PLINK dataset that [`3_Merging_filtering`](../3_Merging_filtering) combines with the LowPass GENCOVE genotypes from [`01_mapping`](../../01_mapping).

**Pipeline engine:** Nextflow **Containerized:** Yes (Seqera Wave) **Estimated runtime:** \[FILL IN\]

------------------------------------------------------------------------

## Pipeline Overview

``` mermaid
%%{init: {'theme': 'neutral'}}%%
flowchart TD
    AXIOM_IN(["Axiom plinkset, canFam3\n(411-dog deposit, or original 804)"])
    IND_KEEP(["411 dogs to keep"])
    CHAIN(["canFam3->canFam4\nliftover chain"])
    REF(["Reference genome"])
    GLIMPSE_REF(["Dog10K phased panel"])
    MAPS(["Genetic maps\n(per chr)"])
    JARS(["conform-gt / beagle jars"])

    LIFTOVER["LIFTOVER_TO_CANFAM4\n(411-sample filter, sex check*,\nliftOver, marker QC*, DA_AFFY)"]
    SPLIT["SPLIT_REF_PANEL_BY_CHR"]
    RECODE["RECODE_AND_RENAME_CHR_VCF\n(+ unused plain .vcf.gz*)"]
    CONFORM["CONFORM_GT"]
    BEAGLE["BEAGLE_IMPUTE"]
    QC["QC_AND_CONVERT_CHR\n(DR2 filter, set IDs,\ncounts*, VCF -> PLINK)"]
    MERGE["MERGE_CHR_PLINK_AXIOM"]

    AXIOM_IN --> LIFTOVER
    IND_KEEP --> LIFTOVER
    CHAIN -.-> LIFTOVER

    LIFTOVER --> RECODE
    REF -.-> RECODE
    RECODE --> CONFORM
    GLIMPSE_REF -.-> SPLIT --> CONFORM
    SPLIT --> BEAGLE
    JARS -.-> CONFORM
    JARS -.-> BEAGLE
    MAPS -.-> BEAGLE
    CONFORM --> BEAGLE --> QC --> MERGE

    OUT[("DA_AFFYimp_ALLCHR .bed/.bim/.fam")]
    MERGE --> OUT
    NEXT(["3_Merging_filtering"])
    OUT --> NEXT
```

`*` marks diagnostic/dead-end sub-steps bundled into that process — reproduced because the original ran them, not because anything downstream consumes their output (see "Notes" below).

------------------------------------------------------------------------

## Input Data

No samplesheet — inputs are direct file/directory params (see [`nextflow_schema.json`](nextflow_schema.json) for the full list).

| Parameter | File | Source |
|-------------------------------|------------------|-----------------------|
| `axiom_bfile_prefix` | `axiom_411_canfam3.{bed,bim,fam}` | Darwin's Ark Axiom array A+B genotypes, canFam3 — deposited as `affy_round3.*` (411 dogs × 1,011,992 variants) at [doi.org/10.17044/scilifelab.33339309](https://doi.org/10.17044/scilifelab.33339309); fetch with `pixi run 00c-fetch-scilifelab-deposit` |
| `ind_to_keep` | `ind_to_keep_round3.txt` | [`data/02_data_processing/ind_to_keep_round3.txt`](../../../data/02_data_processing/ind_to_keep_round3.txt) (411 dogs) |
| `liftover_chain` | `canFam3ToCanFam4.over.chain.gz` | https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/liftOver/canFam3ToCanFam4.over.chain.gz — fetched by `analyses/00_fetch-raw-data/fetch-reference-data.sh` |
| `assembly_ref` | `UU_Cfam_GSD_1.0_ROSY.fa` + `.fai` | Reused from [`01_mapping`](../../01_mapping) |
| `glimpse_ref_panel` | `AutoAndXPAR.Dog10K.phased.bcf` + `.csi` | Reused from [`01_mapping`](../../01_mapping) |
| `chr_start_stop` | `cf4_chr_start_stop.bed` | [`data/02_data_processing/cf4_chr_start_stop.bed`](../../../data/02_data_processing/cf4_chr_start_stop.bed) |
| `chr_rename_map` | `chromosomes.txt` | [`data/02_data_processing/chromosomes.txt`](../../../data/02_data_processing/chromosomes.txt) |
| `genetic_maps_dir` | `mapfiles/canFam4.cM.<chr>.map` | [`data/02_data_processing/mapfiles`](../../../data/02_data_processing/mapfiles) |
| `conform_gt_jar` | `conform-gt.24May16.cee.jar` | https://faculty.washington.edu/browning/conform-gt/conform-gt.24May16.cee.jar — fetch with `pixi run 02a-fetch-axiom-jars` |
| `beagle_jar` | `beagle.22Jul22.46e.jar` | https://faculty.washington.edu/browning/beagle/beagle.22Jul22.46e.jar — fetch with `pixi run 02a-fetch-axiom-jars` |

> The deposited plinkset is the authors' original 804-dog `Affy_merged` already restricted to the 411 dogs in `ind_to_keep_round3.txt` (identical sample IDs) — the output of `LIFTOVER_TO_CANFAM4`'s own first `plink --keep` step — so that step is a no-op on it. Running from the original 804-dog `Affy_merged` (not public) gives the same 411-dog input.

------------------------------------------------------------------------

## Output Data

| Path | Description |
|-----------------------|-------------------------------------------------|
| `plink/DA_AFFYimp_ALLCHR.{bed,bim,fam}` | Merged, imputed Axiom PLINK dataset (canFam4, chr1-38) |
| `vcf/DA_AFFYimp_<chr>.qc.modi.vcf.gz{,.csi}` | Per-chromosome QC'd VCF, pre-merge/pre-PLINK-conversion — consumed by `03_qc`'s Axiom-vs-LowPass concordance study (`bcftools isec`), which needs VCF, not PLINK, format |
| `pipeline_info/versions.yml` | Per-process tool versions |

------------------------------------------------------------------------

## Process Map

| \# | Process | Description |
|----------------|-----------------------|----------------------------------|
| 1 | `LIFTOVER_TO_CANFAM4` | Full canFam3 -\> canFam4 liftover chain: filter to 411 kept dogs, sex check*, build liftover BED, liftOver, derive marker bookkeeping files, filter+flip lifted markers, marker-quality checks*, exclude mismatched-chr markers, update map/chr, re-sort, rename SNP IDs to chr:pos (final `DA_AFFY` plinkset) |
| 2 | `SPLIT_REF_PANEL_BY_CHR` | Materialize one reference-panel chromosome slice |
| 3 | `RECODE_AND_RENAME_CHR_VCF` | Recode one chromosome of `DA_AFFY` to VCF, then rename bare integer chromosome codes to "chrN" (also produces an unused plain-chr `.vcf.gz`\*) |
| 4 | `CONFORM_GT` | Reconcile alleles/strand against the reference panel |
| 5 | `BEAGLE_IMPUTE` | Phase and impute against the Dog10K panel |
| 6 | `QC_AND_CONVERT_CHR` | Filter on imputation quality DR2 \>= 0.8, set variant IDs to CHROM:POS, run diagnostic sample/variant counts\*, convert to PLINK binary |
| 7 | `MERGE_CHR_PLINK_AXIOM` | Merge all chromosomes into the final dataset |

`*` = diagnostic/dead-end sub-step: reproduced because the original ran it, consumed by nothing downstream.

------------------------------------------------------------------------

## Notes

- Each process bundles a group of adjacent original commands that always run together (same scope, no fan-out or other consumer between them) and share a similar resource profile — roughly one process per pipeline stage rather than one per original command. `CONFORM_GT`/`BEAGLE_IMPUTE` are kept as separate processes: `BEAGLE_IMPUTE` needs `process_high` (64 GB, explicit `-Xmx50g`) while the steps around it are `process_low`/`process_medium`, and it is the step most likely to need a memory/retry adjustment. `SPLIT_REF_PANEL_BY_CHR` is also kept separate since the same reference-panel slice is reusable across reruns with different Axiom samples but the same chromosome set.
- `.bim` is used in place of a `.map` file that `01_axiom_liftover.sh`'s liftover-BED-building awk command reads but which is never produced upstream (only `.bed/.bim/.fam` are) — `.bim` carries the same chr/id/cM/pos columns 1-4.
- The sex check, post-liftover marker sanity checks, an intermediate unused plain `.vcf.gz` (`RECODE_AND_RENAME_CHR_VCF`'s `unused_plain_vcf` output), and the sample/variant count checks in `QC_AND_CONVERT_CHR` are diagnostic/dead-end outputs reproduced as leaf outputs even though they feed nothing downstream.
- `LIFTOVER_TO_CANFAM4`'s `A1_notACTG.txt` check reproduces the original awk expression `$5 == !"A/C/T/G"` verbatim; awk evaluates `!"A/C/T/G"` as the boolean negation of a non-empty string (always `0`), so the condition is actually `$5 == 0`, duplicating `A1_is_0.txt` rather than testing for non-ACGT alleles as the filename implies. Preserved as-is.
- Every `--hwe` invocation in this stage writes the `midp` modifier before the p-value threshold (e.g. `--hwe 'midp' 0.000000000000001`), the reverse of plink 1.9's documented `--hwe <p-value> [midp]` order. Preserved as written.
- **Stub test not verified in this pass** (2026-09-18): `pixi run 02b-axiom-imputation -stub` fails schema validation on 5 missing reference/jar paths (`assembly_ref`, `liftover_chain`, `glimpse_ref_panel`, `beagle_jar`, `conform_gt_jar`) — Nextflow's `checkIfExists` runs regardless of `-stub`. Fetchable via `pixi run 00a-fetch-reference-data` + `pixi run 02a-fetch-axiom-jars`, but not fetched here (\~6.75 GB download, \~9 GB on disk after the genome decompresses — skipped for this pass).
- `conform-gt.24May16.cee.jar` and `beagle.22Jul22.46e.jar` (from https://faculty.washington.edu/browning/, fetch with `pixi run 02a-fetch-axiom-jars`) and `canFam3ToCanFam4.over.chain.gz` (from UCSC's `hgdownload.soe.ucsc.edu`, fetched by `analyses/00_fetch-raw-data/fetch-reference-data.sh`) are kept as externally-fetched, version-pinned files rather than conda/bioconda packages: `conform-gt` has no bioconda package, and bioconda's `beagle` currently resolves to a different (2025, `>=5.5_27feb25.75f`) imputation-algorithm release than the `22Jul22.46e` version used here.