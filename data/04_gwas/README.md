# data/04_gwas

Frozen/deposited inputs consumed by [`code/04_gwas`](../../code/04_gwas), plus a few leftover
artifacts from the pre-conversion (pre-Nextflow) workflow. The original, pre-conversion
`README_04.gwas.md` describing the old shell-script pipeline has been moved to
[`archive/legacy_readmes/04_gwas_README_04.gwas.md`](../../archive/legacy_readmes/04_gwas_README_04.gwas.md)
— superseded by each substage's own current `README.md`.

## Live pipeline inputs

| File | Consumed by | Description |
|---|---|---|
| `ALLFAM_DA_MERGED_GENCOVE_AXIOM_QC6.fam` | `3_finemap_susie` (`allfam_fam` param) | Fam file for the full dogCD-GWAS'd cohort |
| `SIZE_DA_MERGED_GENCOVE_AXIOM_QC6.fam` | `3_finemap_susie` (`size_fam` param) | Fam file for the SIZE (body-size) control-trait cohort |
| `STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.fam` | `3_finemap_susie` (`stuck_fam` param) | Fam file for the STUCK (Darwin's Ark phenotype) cohort |
| `SIZE_260603_CSnonfunct.xlsx` | `3_finemap_susie/extract_size_regions.py` | Source spreadsheet SIZE's 9 fine-mapping region windows are extracted from — regenerates the committed `assets/finemap_regions_size.csv` asset; not read at pipeline runtime, only when that asset needs rebuilding |

## Provenance / reference material (not live inputs)

| File | Description |
|---|---|
| `getSIZE_pheno_genofiles.r` | Original pre-conversion source script — `02_data_processing/3_Merging_filtering/bin/build_size_fam.R` reproduces only its SIZE section (see that script's own header comment). Kept for provenance, not run. |
| `NO_PRIORS_finemap_SNPs_in_CS_251122.xlsx` | Frozen output from an earlier (pre-conversion) fine-mapping run — "all SNPs in a credible set," per the archived legacy README. Superseded by the current pipeline's own `finemap_results` output; kept for reference. |

## Reader-facing summary tables (not a live pipeline input, kept per KatarinaTe)

| File | Description |
|---|---|
| `ALLFAM_sex_age_phenos.txt` | Per-dog table (2,585 rows): dog ID, age, sex, F1/F2/F3 factor scores + SEs, all 14 item responses, for the full ALLFAM cohort |
| `ALLFAM_sex_age_phenos_batch.txt` | Same table, with an added `batch` column |
| `sex_batchCCDF1_DA_MERGED_GENCOVE_AXIOM_QC6.covar` | Plink-style covariate file (2,429 rows) restricted to the F1-specific QC6 cohort |

Not read by any current script — matches variable/file names that are **commented out** in the
original source script (`03_create_files_per_factor_mlma.R`'s `#write.table(...)` calls), so the
live pipeline doesn't regenerate these either. Kept intentionally: these are intermediate
summaries provided to readers, not orphaned scratch output.

These are exploratory covariate tables from the F1 batch-covariate calibration check (does adding
sequencing/genotyping `batch` as a GWAS covariate change calibration?) — the underlying data
behind the manuscript's "batch-covariate-in-GWAS comparison" text (flagged in
`REVIEW_AND_SUGGESTIONS.md`'s open items as needing to move from Results to Supp Info). The actual
analysis — GCTA run twice with different covariates, then compared — lives in
[`code/04_gwas/1_mlma-loco/`](../../code/04_gwas/1_mlma-loco) (`MLMA_PER_FACTOR.sh`,
`compare_batch_gwas_calibration.R`); see that stage's own README for exact run instructions.
