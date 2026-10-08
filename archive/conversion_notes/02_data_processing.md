# Conversion notes — 02_data_processing

Decisions and discrepancies noted while converting this stage's original scripts into the
Nextflow pipelines under `code/02_data_processing/`. Moved out of the substage `README.md` files
and `bin/*` script headers so those files document only the science/pipeline itself.

## code/02_data_processing/1_Survey_data/README.md

### Design Notes

#### Every diagnostic and plot is reproduced as a real output file

The original script was run interactively — a bare `plot()`/`ggplot()` call renders to an
interactive graphics device, and a bare top-level expression (e.g. `fa.kmo`, `head(x)`) auto-prints
to the console. Neither happens automatically in a non-interactive container run, so every such
call in `bin/*.R` is wrapped: plots in `pdf(...)`/`dev.off()`, and text diagnostics inside a
`sink(log_file, split = TRUE)` block spanning the whole computational section (so nothing has to be
individually enumerated and risk being missed). This mirrors the "reproduce everything the author
ran, not just what feeds the final output" principle already applied in `2_Axiom_imputation` — here
it applies to essentially the entire script, since almost none of `RUN_EFA`'s output and much of
each factor's diagnostic output feeds anything downstream.

#### Per-factor items are derived by column, not re-derived from raw survey data

The original's F1/F2/F3 sections each rebuild their item table from the raw long-format answers
data (not in this repo) and re-apply the same "switch"/"rescale" recoding already baked into
`response_df_dogCDitems.txt`. Checked column-by-column: each factor's recoding only ever touches
that factor's own items, independent of what else is in the table — so selecting columns from the
already-recoded checkpoint produces byte-identical values to re-deriving from raw data.
`RUN_IRT_FACTOR` does the former.

#### Real per-factor differences preserved, not normalized into a uniform loop

`RUN_IRT_FACTOR` is one config-driven process (avoiding the DSL2 "process already used" issue)
rather than three near-identical processes, but the per-factor differences in the original are
genuine and are passed through as config, not smoothed away:

- `M2()`'s `type` argument: F1 and F3 use `"M2*"`, F2 uses `"C2"`.
- Which `fscores()` method feeds the saved CSV: F1 and F3 use MAP, F2 uses EAP (matching each
  block's own "combine with dog ID — using MAP/EAP" comment). Both EAP and MAP are always computed
  and plotted for every factor regardless — only the one joined into the final output differs.
- F1 has one extra diagnostic scatter plot (`item7` vs. factor SE) that F2/F3 don't have.
- `tracePlot`/`testInfoPlot`/`itemInfoPlot`'s `theta_range` and `scale_color_brewer`'s palette are
  distinct per-factor values, not shared cosmetic defaults — passed through as config.

#### A few small per-factor asymmetries are not reproduced

A handful of tiny diagnostic-log print differences between F1/F2/F3 (an extra `head()`/`nrow()`
here, a bare filter print there, `head()` vs `tail()` on the saved data.frame) and the EAP/MAP
plot's exact per-factor `title()`/`xlim`/point color are not reproduced — none of them affect the
actual computed scores or the saved `CCD3F.txt`, only incidental console-inspection ordering and
plot cosmetics. The EAP/MAP and NA-count plots use one consistent generic form
(`title = paste0(factor_id, "_CCD3F")`) across all three factors instead of chasing each factor's
exact original styling.

### Reproducibility Caveats

#### `nFactors::parallel()` is unseeded

The random parallel-analysis simulation (`rep=100`) that produces the `nScree` diagnostic plot has
no set seed, so it isn't reproducible run-to-run. This only affects that one diagnostic plot,
though — the actual number of factors used (`nFac = 3`) is a fixed value already hardcoded in the
original script, and the EFA fit itself (`fa(..., nfactors = 3, fm = "pa")`) is a deterministic
iterative fit given fixed inputs, independent of the random simulation. Not seeded here either,
matching the original exactly rather than "fixing" it to something more reproducible than what
actually ran.

#### Why 3 factors, not 4

Verified against real data (`fa(..., nfactors = 4)`): the 4th factor's SS loadings = 0.669 (below
1, the standard retention threshold). With `nfactors = 3`, all three factors have SS loadings above
1 (2.126, 1.357, 1.063).

#### `mirt()`'s IRT fitting

Flagged as a reproducibility concern in the project audit alongside the parallel-analysis issue
above. Not independently verified here — R isn't available in this development environment to
check empirically whether the default EM estimation used is actually stochastic or was flagged out
of caution. Documented rather than asserted either way.

#### Hardcoded `item.weights`

Each factor's `item.weights` vector (e.g. F1: `c(0.81,4.86,25.77,1.44)`) is a literal in the
original script. It looks copy-pasted from an earlier `mirt()` run's printed `params$items`
discrimination parameters rather than recomputed each run — kept as a hardcoded value here
(`FACTOR_CONFIGS` in the main workflow) rather than recomputed from `params$items` at runtime,
matching how the original actually used it.

### Resolved Issues from Original Implementation

#### Redundant `nfactors()` call removed from `RUN_EFA` (a scoped exception, not general policy)

The original calls `nfactors(response_df3, n=10)` on the not-yet-NA-omitted data, then immediately
overwrites the result with `nfactors(response_df4, n=10)` on the NA-omitted version two lines
later — the first call's result is never used or displayed, and it operates on data that's about
to be superseded for every other diagnostic in the script. Removed in `run_efa.R` as a genuine
pointless duplicate, unlike the many other preserved-but-unused diagnostics elsewhere in this
pipeline (e.g. `2_Axiom_imputation`'s `CHECK_SEX`, or this script's own double `summary(fitGraded)`
calls in `RUN_IRT_FACTOR`) — this is a one-off cleanup specific to this exact redundant pair, not a
change to the general "reproduce everything the author ran" policy.

#### `M2()` always errors — reproduced as a caught, logged error, not a script crash

Every one of F1/F2/F3's `M2()` calls carries the original's own comment `#Error: M2() statistic
cannot be calculated due to too few degrees of freedom` — the original author documented hitting
this exact error every time, for every factor; it's a structural fact about this model, not
data-dependent noise. Run interactively, an error on one line just prints and the author moved on
to the next line by hand; a batch `Rscript` run halts the entire script on an uncaught error
instead. `run_irt_factor.R` wraps the call in `tryCatch()` so the error is logged (matching what
the interactive session actually showed) and execution continues to `itemfit()` and beyond —
verified against real data, confirmed to reliably reproduce this exact error for all three
factors, and confirmed the `tryCatch` correctly lets the rest of the script run to completion.

#### Missing `ggplot2` import

`run_irt_factor.R` composes `ggmirt`'s `tracePlot()`/`itemInfoPlot()` output with `+ labs(...)` and
`+ scale_color_brewer(...)` (`ggplot2` functions), matching the original script's own code exactly
— but the original relies on `ggplot2` already being attached globally via its `library(tidyverse)`
call at the very top of the file. Since this pipeline only loads what's used past the
`response_df_dogCDitems.txt` checkpoint (see `env/02_data_processing/README.md`), `ggplot2` needed
an explicit `library()` call here that the original never needed to state directly. Caught by
running the real (non-stub) pipeline against real data, not by inspection.

## code/02_data_processing/2_Axiom_imputation/README.md

### Process Granularity

Each process above bundles a group of adjacent original commands that always run together (same
scope, no fan-out or other consumer between them) and share a similar resource profile —
roughly one process per pipeline stage rather than one per original command. Steps were **not**
merged across a fan-out boundary, and `CONFORM_GT`/`BEAGLE_IMPUTE` were kept separate from each
other and from their neighbors despite running back-to-back per chromosome, because:

- `BEAGLE_IMPUTE` needs `process_high` (64 GB, explicit `-Xmx50g`) while every step around it is
  `process_low`/`process_medium` — merging would mean requesting the larger allocation for the
  whole merged task's duration, including the cheap parts.
- `BEAGLE_IMPUTE` is the step most likely to need a memory/retry adjustment on real data (the
  original author's own comments note having to bump its heap size mid-project) — keeping it its
  own process means a retry there doesn't also re-run `CONFORM_GT`.
- `SPLIT_REF_PANEL_BY_CHR` is an independently reusable/cacheable concern (the same reference-panel
  slice would be valid across reruns with different Axiom samples but the same chromosome set), so
  it's kept out of the steps that consume it.

Merging elsewhere (the 11-step liftover chain into `LIFTOVER_TO_CANFAM4`; `RECODE_CHR_VCF` +
`RENAME_CHR_VCF` into `RECODE_AND_RENAME_CHR_VCF`; `QC_FILTER_DR2` + `SET_VARIANT_IDS` +
`COUNT_IMPUTED_VARIANTS` + `VCF_TO_PLINK` into `QC_AND_CONVERT_CHR`) costs nothing in lost
parallelism or `-resume` granularity that would actually matter, and saves a proportional number of
container starts per run (21 processes → 7).

### Resolved Issues from Original Implementation

#### `.map` file referenced but never produced

`01_axiom_liftover.sh`'s liftover-BED-building `awk` command reads `affy_round3.map`, but the
preceding `plink --make-bed` only ever produces `.bed/.bim/.fam` — no `.map` is generated anywhere
upstream. `.bim` carries the same chr/id/cM/pos columns 1-4 that the awk command reads, so
`LIFTOVER_TO_CANFAM4` uses `.bim` in its place.

#### Diagnostic and dead-end steps: reproduced, not omitted

Several original commands produce output that nothing downstream ever reads: the `--split-x
--check-sex` sex check and the post-liftover marker sanity checks on `affy_cf3.bim` (both bundled
into `LIFTOVER_TO_CANFAM4`), the intermediate `bcftools view DA_AFFY_$chrN.vcf -Oz -o
DA_AFFY_$chrN.vcf.gz` created right before the chromosome-rename step re-reads the plain `.vcf`
instead (`RECODE_AND_RENAME_CHR_VCF`'s `unused_plain_vcf` output), and the sample/variant count
checks at the end of `02_axiom_impute.sh` (bundled into `QC_AND_CONVERT_CHR`). All four are
reproduced as leaf outputs/stdout-only commands even though they feed nothing else — the goal is to
redo what the original author ran, not just the subset of commands that happen to feed the final
output.

#### `$5 == !"A/C/T/G"` in the non-ACGT marker check

`LIFTOVER_TO_CANFAM4`'s `A1_notACTG.txt` check reproduces the original awk expression
`$5 == !"A/C/T/G"` verbatim. awk evaluates `!"A/C/T/G"` as the boolean negation of a non-empty
string, which is always `0` — so the condition is actually `$5 == 0`, duplicating `A1_is_0.txt`
rather than testing for non-ACGT alleles as the filename implies. Preserved as-is, not corrected.

#### `--hwe 'midp' <p>` argument order

Every `--hwe` invocation in the `02_data_processing` scripts writes the `midp` modifier before the
p-value threshold (e.g. `--hwe 'midp' 0.000000000000001`), the reverse of plink 1.9's documented
`--hwe <p-value> [midp]` order. Reproduced literally: every archived run reports the expected
post-filter variant counts immediately after, confirming the command ran and filtered as intended
regardless of argument order — this is preserved, not corrected.

#### Liftover chain and imputation jar provenance

Neither `canFam3ToCanFam4.over.chain.gz` nor the `conform-gt`/`beagle` jars were documented
anywhere in this repo with a source. Tracked down and pinned to the exact releases named in
`02_axiom_impute.sh`'s comments: `conform-gt.24May16.cee.jar` and `beagle.22Jul22.46e.jar` from
https://faculty.washington.edu/browning/ (fetch with `pixi run 02a-fetch-axiom-jars`), and
`canFam3ToCanFam4.over.chain.gz` (dated 2020-05-12, matching the script's comment) from UCSC's
`hgdownload.soe.ucsc.edu` (fetched by `analyses/00_fetch-raw-data/fetch-reference-data.sh`). Both
jars are kept as externally-fetched, version-pinned files rather than conda/bioconda packages:
`conform-gt` has no bioconda package at all, and while bioconda does package `beagle`, it currently
resolves to a 2025 build (`>=5.5_27feb25.75f`) — a different imputation-algorithm version than the
`22Jul22.46e` release actually used, which would change results rather than merely repackage them.

## code/02_data_processing/4_Plotting/README.md

### Design Notes

#### Every diagnostic and plot is reproduced as a real output file

Like `1_Survey_data` and `2_Axiom_imputation`, the original scripts were run interactively — bare
`plot()`/`ggplot()` calls render to an interactive graphics device, and bare top-level expressions
(`head(x)`, `table(x)`, `nrow(x)`) auto-print to the console. Every such call in `bin/*.R` is
wrapped: plots in `png(...)`/`dev.off()`, text diagnostics inside a `sink(log_file, split = TRUE)`
block.

#### Genericized via per-plot config tables, not near-duplicated ggplot blocks

`plot_factor_distr.R` and `plot_item_factor_hist.R` build 3 and 17 plots respectively from small
config tables (one row per plot) rather than repeating near-identical ggplot code per plot. Every
functional/styling difference the original actually has — title, axis label, color, x-axis
breaks/limits, title/axis-title margins, axis hjust, whether the y-axis is blanked — is preserved
exactly per plot, per direct user request for exact-margin fidelity (2026-09-01), including one
genuine source oddity kept as-is: item152 has no `scale_x_continuous(breaks=1:6)` call unlike its
6 CCDF3 siblings.

#### Hardcoded values baked into the visible plots are computed from real data instead

Several values in the original are hardcoded numbers baked into the *visible* plot output, not
just comments — PCA's axis-label variance percentages, every `N=` annotation, and the shared axis
ranges for the continuous factor-score histograms — all clearly correct for whatever historical
run the original author had on hand, but not guaranteed to match a fresh run's real counts/ranges.
Per direct user instruction (2026-09-01), these are computed from the real data here instead: N
via `nrow()`/`sum(!is.na())`, PCA variance-% from the real eigenval file, and axis ranges from the
real combined data range (with a small pad) rather than the original's specific hardcoded numbers.
The 14 items' fixed 1-6 Likert response-scale breaks are the one exception kept literal — that's a
survey-instrument property, not a stale historical count.

#### Container reused, not rebuilt

`4_Plotting` needs exactly the same R/dplyr/ggplot2/patchwork/ggtext stack `3_Merging_filtering`
already has frozen — see [`env/02_data_processing/README.md`](../../../env/02_data_processing/README.md).

### Reproducibility Caveats

#### `dog_breeds_csv` — external file, now resolved

The original reads breed metadata from a literal placeholder path
(`/path/to/DarwinsArk_20220715_dogs.csv`) — a real path scrubbed before the script was committed.
The real file has since been added to the repo
(`data/02_data_processing/DarwinsArk_20220715_dogs.csv`) and is wired up as the required
`dog_breeds_csv` param, validated against it directly (columns confirmed: `dog`, `purebred`,
`breed1_inputted`).

#### PCA's uncolored first plot is dropped, not reproduced

The original's first PCA plot has no color at all (`col=col.vector` is commented out in its own
`plot()` call), yet its legend still lists breed colors ("other breed"/"Labrador Retriever" in
red/blue) — a real, pre-existing inconsistency in the original script. Per direct user
instruction (2026-09-01), this plot isn't reproduced at all; only the second, fully colored/
shaped-by-breed plot (`pca_by_breed.png`) is kept.

#### item147's `limits = c(0.5, 6, 5)` fixed to `c(0.5, 6.5)`

A 3-element vector where `scale_x_continuous(limits=...)` expects 2 (verified harmless on its
own — ggplot silently uses only the first two elements). Fixed here per direct user observation:
this codebase has other Swedish-locale tells (e.g. `SD2VÄRDE` in the original `plot_factor_distr.R`),
making "6,5" typed as a European decimal comma — silently split by R into two vector elements
inside `c(...)` instead of one `6.5` — a far more likely explanation than a deliberately
open-ended upper bound.

## code/02_data_processing/1_Survey_data/process_definitions.nf

Original file header:

```
// ============================================================================
// process_definitions.nf — 1. Survey data (EFA/IRT)
//
// Source: 02_data_processing/1_Survey_data/Darwin_dog-CD_3-factorsolution.R
//
// Every diagnostic print and plot in the original script is reproduced as a real output file
// here — sink()-captured logs and png()-captured plots — rather than left to whatever the
// interactive R console happened to show, matching the same "reproduce everything, not just what
// feeds the final output" principle used in 2_Axiom_imputation. See bin/*.R for exactly which
// original lines each output corresponds to. Plots are PNG, not the original's PDF, per direct
// user request (2026-08-28) — RUN_IRT_FACTOR's na_count_plots (originally 3 pages in one PDF)
// splits into 3 separate PNGs accordingly; see its header comment in bin/run_irt_factor.R.
// ============================================================================
```

Original `RUN_EFA` comment:

```
    // EFA diagnostics (nfactors/parallel analysis/nScree, KMO, Bartlett, determinant, cronbach's
    // alpha) and fit. Source: lines 119-236. None of this feeds the F1/F2/F3 IRT models below —
    // the item groupings they use are a fixed, pre-decided grouping, not derived programmatically
    // from this fit.
```

Original `RUN_IRT_FACTOR` comment:

```
    // Per-factor IRT fit and factor scores. Source: F1 lines 239-356, F2 lines 357-469,
    // F3 lines 470-577 — genericized via arguments, not three near-duplicate processes. See
    // bin/run_irt_factor.R's header for the real per-factor differences preserved (M2() `type`,
    // which fscores() method feeds the saved CSV, F1's one extra diagnostic plot, and the
    // tracePlot/testInfoPlot/itemInfoPlot theta ranges + color palette, none of which are uniform
    // across factors despite looking like cosmetic plot styling).
```

Original `PLOT_EXTRA_ITEM_SE` comment:

```
    // F1's one extra diagnostic scatter plot (an item's response vs. factor SE), that F2/F3
    // don't have. Source: lines 347-348. Split out from RUN_IRT_FACTOR and invoked only for F1
    // (see the main workflow) rather than forcing RUN_IRT_FACTOR to declare a conditionally-
    // produced output — Nextflow output paths require a minimum arity of 1, so "produced only
    // for some factors" can't be expressed as an optional output on a process all three factors
    // share. Reads the factor's already-saved CCD3F.txt rather than duplicating any computation.
```

Original `PLOT_FACTOR_DISTRIBUTIONS` comment:

```
    // Final per-factor score distribution plots. Source: lines 579-615. Needs all three factors'
    // CCD3F files together — a barrier after the per-factor fan-out.
```

## code/02_data_processing/2_Axiom_imputation/process_definitions.nf

Original `LIFTOVER_TO_CANFAM4` comment:

```
    // Full canFam3 -> canFam4 liftover chain for the Axiom array genotypes.
    // Source: 02_data_processing/2_Axiom_imputation/01_axiom_liftover.sh
    //
    // Merged from 11 originally-separate steps (filter to 411 samples, sex-check diagnostic,
    // build liftover BED, liftOver, derive marker bookkeeping files, filter+flip lifted markers,
    // marker-quality diagnostics, exclude mismatched-chr markers, update map/chr, re-sort
    // positions, rename SNP IDs): a strictly linear, cohort-level, single-plinkset chain with no
    // fan-out and no other consumer of any intermediate file, so merging it costs nothing in lost
    // parallelism or resume granularity, and saves 11 container starts per run.
    //
    // Resolved issues, still applicable:
    // - `.map` file referenced but never produced: the awk command reads affy_round3.map, but
    //   `plink --make-bed` only produces `.bed/.bim/.fam`. `.bim`'s chr/id/cM/pos columns 1-4
    //   match a `.map` file's, so `.bim` is used in its place.
    // - `$5 == !"A/C/T/G"` in the marker-quality check: preserved verbatim. awk evaluates this as
    //   `$5 == 0`, duplicating the A1_is_0 check rather than testing for non-ACGT alleles as its
    //   filename implies. Not corrected.
    // - `--hwe 'midp' <p>`: argument order preserved as written (see the pipeline README).
    // - The sex check and marker-quality checks are diagnostics whose output nothing downstream
    //   reads — reproduced anyway, since the goal is to redo what the original author ran, not
    //   just the subset of commands that feed the final output.
```

Original `QC_AND_CONVERT_CHR` comment excerpt:

```
    // The diagnostic counts read a mix of pre-QC (the raw beagle output) and post-QC (qc.modi)
    // files, matching the original exactly — see the equivalence walkthrough. Reproduced as
    // stdout only, captured in this task's own log.
```

## code/02_data_processing/4_Plotting/process_definitions.nf

Original file header:

```
// ============================================================================
// process_definitions.nf — 4. Plotting
//
// Source: 02_data_processing/4_Plotting/{plot_PCA.R, plot_factor_distr.R,
//          plot_item&factor_hist.R}
//
// Reuses 3_Merging_filtering's container (community.wave.seqera.io/library/
// merging_filtering_tools) — same R/dplyr/ggplot2/patchwork/ggtext stack, no new build needed.
//
// Every diagnostic print in the original scripts is reproduced as a real sink()-captured log —
// see bin/*.R for exactly which original lines each output corresponds to. Several hardcoded
// values baked into the *visible* plots in the original (PCA variance-% axis labels, N=
// annotations, axis ranges for continuous scores) are computed from the real data here instead,
// per direct user instruction (2026-09-01) — a plot's own numbers should never lie about the data
// actually plotted. See each bin/*.R header for exactly what's computed vs. kept literal.
// ============================================================================
```

Original `PLOT_PCA` comment:

```
    // PCA scatter plot colored/shaped by breed. Source: plot_PCA.R (full script). Breed
    // metadata (dog_breeds_csv, data/02_data_processing/DarwinsArk_20220715_dogs.csv) — see
    // bin/plot_pca.R's header for the variance-% axis-label fix and why the original's uncolored
    // first plot (whose legend didn't match its own data) is dropped rather than reproduced.
```

Original `PLOT_FACTOR_DISTR` comment:

```
    // F1/F2/F3 factor-score distribution histograms, stacked vertically. Source:
    // plot_factor_distr.R (full script). See bin/plot_factor_distr.R's header for the
    // per-factor config (only F3's block styles axis.title.x in the original) and the
    // dynamically-computed N=/shared-axis-range fix.
```

Original `PLOT_ITEM_FACTOR_HIST` comment excerpt:

```
    // Factor-score histogram row (F1+F2+F3) and three per-item response histogram rows (CCDF1's
    // 4 items, CCDF2's 3 items, CCDF3's 7 items). Source: plot_item&factor_hist.R (full script).
    // See bin/plot_item_factor_hist.R's header for the per-plot config tables (title, color,
    // x-axis breaks/margins — reproduced exactly per plot, including item152's missing
    // scale_x_continuous() call, a genuine source oddity kept as-is) and the item147
    // `c(0.5,6,5)` -> `c(0.5,6.5)` fix (a likely European-locale decimal-comma typo, per direct
    // user observation — this codebase has other Swedish-locale tells, e.g. SD2VÄRDE).
```

## code/02_data_processing/3_Merging_filtering/process_definitions.nf

Original `FLIP_SCAN_CONVERGENCE` comment:

```
process FLIP_SCAN_CONVERGENCE {
    // Iterative strand-flip detection and exclusion, converging when a round's flip-scan finds
    // nothing left to exclude. The archived run took exactly 2 rounds of real exclusion (7163,
    // then 264 variants) before a 3rd scan found zero — this loop is written to converge
    // generically rather than assuming exactly that count, since it's a purely file-based
    // criterion (an empty exclude list) with no human judgment involved, unlike the KING-based
    // duplicate/relatedness exclusion later in this pipeline.
    //
    // Single process with an embedded bash loop, not one Nextflow process per round: DSL2
    // processes can't be invoked in a loop, and there's no fan-out here to lose by keeping it as
    // one task.
    //
    // Two things resolved, not reproduced literally:
    // - Round 1's `--make-pheno DA_IMP_GENCOVE.fam '*'` references a file never created anywhere
    //   in these scripts (only `DA_IMP_GENCOVE_QC.fam` exists, which rounds 2+ correctly use) —
    //   standardized to `DA_IMP_GENCOVE_QC.fam` for every round.
    // - Round 1's `--freq` diagnostic checks on the flipped markers (comparing MAF between the two
    //   platforms) are not reproduced — round 2 has no equivalent, and this is diagnostic-only,
    //   feeding nothing downstream.
    //
    // One thing added that the original never needed: `.flipscan` has a header row, and
    // `awk '$10 != "NA"'` alone would keep that header line forever (its 10th field is a column
    // label, never literally "NA"), so the derived exclude list would never be truly empty and
    // this loop would never converge. The original never hit this because it never relied on
    // emptiness as a stopping condition — it ran a fixed 3 rounds and stopped by hand once the
    // printed summary said 0 hits. Automatic convergence needs `NR>1` to skip the header;
    // verified with a mocked plink/awk harness that the loop hangs without it and converges
    // correctly with it.
```

Original `APPLY_EXCLUDE_LOWPASS12` comment:

```
process APPLY_EXCLUDE_LOWPASS12 {
    // Not in the original script. 01_mapping's own README documents 12 samples as "excluded in
    // the end (02_data_processing/3_Merging_filtering/01_gencove_axiom_merging_filtering.sh)"
    // due to no Darwin's Ark dogID and no phenotypes — but that script never actually removes
    // them: their missing depth value makes build_exclude_list_round1.R's DEPTH_filter2 default
    // to 1 ("keep"), same as every other axiom sample with no depth data, and none of the three
    // exclude-list files name them either. Confirmed they ride through QC4C untouched (verified
    // present in modi_DA_MERGED_GENCOVE_AXIOM_QC4C_forDogIDchange.fam and in the frozen data6
    // checkpoint as all-NA rows). Added per direct user instruction after that investigation.
    //
    // Removed here from QC4C's own native fam (still keyed by raw sampleID at this point — the
    // dogID relabel hasn't happened yet) via plink's ordinary ID-based --remove, the same
    // mechanism as APPLY_EXCLUDE_ROUND1/2 — not by trimming the dogID-relabel crosswalk fam
    // directly, which would desync its row order from QC4C's .bed/.bim before
    // RELABEL_AND_EXCLUDE_ROUND3's positional --fam substitution (plink matches --fam rows to
    // .bed/.bim rows by position, not by ID).
    //
    // exclude_list_lowpass12 (exclude_12_LowPass_samples.txt) has a header row and a single ID
    // column; reformatted inline to plink's headerless two-column FID/IID --remove format
    // (FID==IID for these samples, same convention list_exclude.txt already uses).
    //
    // Scoped deliberately to only these 12: a broader check (data/02_data_processing/
    // DA_MERGED_GENCOVE_AXIOM_QC3modi.fam rows where sampleID==dogID, cross-verified against
    // absence from the covariates/depth files) found 7 more samples with an identical profile
    // (2 Axiom "-a" secondary-run IDs, 5 GENCOVE SRR accessions) that are NOT named in any
    // provided exclude list. Flagged for the user, not folded in here.
```

Original `RELABEL_AND_EXCLUDE_ROUND3` comment excerpt:

```
    // Inputs are now the LowPass-12-filtered QC4C and crosswalk fam (see APPLY_EXCLUDE_LOWPASS12/
    // FILTER_DOGID_RELABEL_FAM above) rather than QC4C/the crosswalk fam directly — a deliberate
    // deviation from the original, whose own QC5 dog count (3328) will no longer match: expect
    // 3328 - 12 = 3316 dogs here instead (none of the 12 overlap the dog-3094 duplicate pair
    // this process also removes).
```

Original `BUILD_SIZE_COVARIATES` comment excerpt:

```
    // filenames hardcode a "CCD" prefix, but SIZE's real files don't use one (same reasoning as
    // FILTER_SIZE_QC6). Reuses bin/build_factor_covariates.R unchanged — that script takes its
    // output prefix as a plain argument, already fully generic. Per direct user confirmation
    // (2026-09-01): SIZE's real historical mlma-loco run used age/sex covariates built exactly
    // like the factors.
```

Original Stage E file-header block:

```
// ============================================================================
// Stage E: per-factor phenotype files (F1/F2/F3)
// Source: 02_data_processing/3_Merging_filtering/03_create_files_per_factor_mlma.R
//
// `data6_2024-10-14.txt` is consumed here as a frozen input rather than recomputed from data1 +
// F1/F2/F3_CCD3F.txt (the join the original script itself performs immediately before reading
// this file back in, overwriting the freshly-computed result — the same
// "write-commented-then-read-back" idiom resolved elsewhere in this pipeline by recomputing
// fresh). Not resolved that way here: the original's `age = data5$age.x, sex = data5$sex_numeric.x`
// column references imply a name collision between data1 and the per-factor CCD3F tables that
// doesn't reproduce with the CCD3F files this pipeline's own 1_Survey_data stage produces (which
// carry no `age`/`sex_numeric` columns to collide with) — recomputing would silently drop those
// columns (or error on first downstream use), and there is no committed historical F1/F2/F3
// CCD3F file to check the original collision against. Flagged for the user rather than resolved
// unilaterally; using the frozen, already-verified-nonempty checkpoint sidesteps the ambiguity.
//
// data6 also needs restricting to QC5's own dogID set/order before BUILD_FACTOR_FAM/
// BUILD_ITEM_FAM can use it — see FILTER_DATA6_TO_QC5's header for why (a consequence of the
// LowPass-12 fix in Stage C). BUILD_ITEM_GRAB (Stage F) deliberately keeps using the
// *unfiltered* data6, since its join there is ID-based, not a positional --fam substitution, and
// that's what the original literally does.
// ============================================================================
```

Original `FILTER_DATA6_TO_QC5` comment:

```
process FILTER_DATA6_TO_QC5 {
    // Restricts + reorders the frozen data6 checkpoint to exactly QC5's own dogID set, in QC5's
    // own row order. See bin/filter_data6_to_qc5.R for the full rationale: QC5 is 3316 rows
    // after Stage C's LowPass-12 fix, data6 is still 3328 (frozen, untouched by that fix), and
    // BUILD_FACTOR_FAM/BUILD_ITEM_FAM's output feeds a strictly positional plink --fam
    // substitution downstream (FILTER_FACTOR_QC6/FILTER_ITEM_QC6) that requires exact row-count
    // and row-order agreement with QC5's .bed/.bim.
```

Original SIZE file-header banner:

```
// ============================================================================
// SIZE — control-phenotype fam file, run via mlma alongside F1/F2/F3 (not part of the
// original 01/02/03/04 shell/R sources for this pipeline — see data/04_gwas/
// getSIZE_pheno_genofiles.r). Added per direct user instruction (2026-09-01). That file's own
// second section (STUCK, a different Darwin's Ark item) is deliberately not reproduced here.
// ============================================================================
```

## code/02_data_processing/3_Merging_filtering/3_merging_filtering.nf

Original file header:

```
// ============================================================================
// main.nf — 3. Merging & filtering workflow
//
// Merges the LowPass GENCOVE (01_mapping) and Axiom (2_Axiom_imputation) genotype datasets,
// filters samples and variants, and produces per-phenotype (3 factors + 14 items) GWAS-ready
// PLINK datasets.
//
// Source: 02_data_processing/3_Merging_filtering/{01_gencove_axiom_merging_filtering.sh,
//          02_create_ALLFAM_check_stats.sh, 03_create_files_per_factor_mlma.R,
//          04_create_files_per_item_polmm.R}
//
// Built incrementally, stage by stage — implements Stage A (cohort QC + merge), Stage B
// (flip-scan convergence), Stage C (sample-level QC: duplicates, depth, KING relatedness,
// dogID relabeling), Stage D (ALLFAM stats, PCA, GRM), Stage E (per-factor phenotype files),
// Stage F (per-item phenotype files), and Stage G (cross-phenotype dog/SNP-count check).
//
// One deliberate deviation from the original within Stage C: APPLY_EXCLUDE_LOWPASS12 /
// FILTER_DOGID_RELABEL_FAM remove 12 samples that 01_mapping's own README documents as excluded
// here, but which 01_gencove_axiom_merging_filtering.sh never actually removes. Added per direct
// user instruction — see APPLY_EXCLUDE_LOWPASS12's header comment for the full investigation.
//
// A second, consequent deviation at the Stage C/E boundary: that same fix shrinks QC5 to 3316
// dogs, while the frozen data6 checkpoint (untouched by it) still has 3328 — so data6 is
// restricted + reordered to QC5's own dogID set/order (FILTER_DATA6_TO_QC5) before
// BUILD_FACTOR_FAM/BUILD_ITEM_FAM use it. See FILTER_DATA6_TO_QC5's header comment.
// ============================================================================
```

## code/02_data_processing/1_Survey_data/bin/run_efa.R

Original header:

```
#!/usr/bin/env Rscript
# run_efa.R — EFA diagnostics and fit for the 14 dogCD survey items
# Source: 02_data_processing/1_Survey_data/Darwin_dog-CD_3-factorsolution.R, lines 119-236
#
# All diagnostic prints in this script (nfactors/KMO/Bartlett/det/cronbach's alpha/loadings) feed
# nothing downstream in the original — the factor solution actually used (F1/F2/F3's item lists)
# is a fixed, pre-decided grouping, not something this script derives programmatically. Captured
# anyway: the goal is to redo everything the original author ran, not just what feeds the final
# output.
```

Original inline comment (before the NA-omit step):

```
# remove NAs to continue
#
# The original also runs `nfactors(response_df3, n=10)` here first, on the not-yet-NA-omitted
# data, immediately before this line overwrites it with the NA-omitted version — that call's
# result is never used or displayed, and every other diagnostic in this script operates on the
# NA-omitted response_df4 (or mat_cor derived from it), never on response_df3 directly again.
# Removed here as a pointless duplicate computation on soon-to-be-superseded data — a
# process-specific cleanup, not a general policy of dropping preserved-but-unused originals.
```

## code/02_data_processing/1_Survey_data/bin/run_irt_factor.R

Original header:

```
#!/usr/bin/env Rscript
# run_irt_factor.R — per-factor IRT fit and factor scores
# Source: 02_data_processing/1_Survey_data/Darwin_dog-CD_3-factorsolution.R
#   F1: lines 239-356, F2: lines 357-469, F3: lines 470-577
#
# Genericized across the three factors via arguments rather than three near-duplicate scripts.
# Real per-factor differences preserved, not normalized away:
#   - M2()'s `type` argument: F1 and F3 use "M2*", F2 uses "C2".
#   - Which fscores() method feeds the saved output: F1 and F3 use MAP (tabscores3), F2 uses EAP
#     (tabscores1) — both EAP and MAP are always computed and plotted for every factor; only the
#     one that gets joined into the saved CSV differs.
#   - F1 has one extra diagnostic scatter plot (an item's response vs. factor SE) that F2/F3
#     don't — handled by a separate process (PLOT_EXTRA_ITEM_SE), invoked only for F1, reading
#     this script's saved CSV output rather than duplicating any computation here.
#   - Plot ranges/palette genuinely differ per factor, not just item count: tracePlot's
#     theta_range, testInfoPlot's theta_range, itemInfoPlot's (non-faceted) theta_range and
#     scale_color_brewer palette are all distinct per-factor values in the original, not cosmetic
#     styling — passed through as config rather than hardcoded.
#
# Not reproduced (see the pipeline README): a handful of small per-factor asymmetries in
# diagnostic-log print ordering/duplication (e.g. an extra head()/nrow() here, a bare filter
# print there) and in the EAP/MAP plot's exact title/xlim/color — none of which affect the actual
# computed scores or the saved CCD3F.txt. The EAP/MAP and NA-count plots use one consistent
# generic form across all three factors instead.
#
# response_df2 is subsetted by column from the already-recoded, already-committed
# response_df_dogCDitems.txt checkpoint rather than re-derived from the raw long-format answers
# table (which isn't in this repo). This is a verified equivalence, not an approximation: each
# factor's original "switch"/"rescale" recoding only ever touches that factor's own items, so
# selecting columns from the already-recoded 14-item table gives byte-identical values to
# re-deriving from raw data (see the pipeline README).
#
# Column selection by name (`item_cols`) replaces the original's per-factor hardcoded positional
# indices (e.g. F1's `x[,2:5]`, F2's `x[,2:4]`, F3's `x[,2:8]`) — those ranges exist only because
# each factor has a different item count; selecting by name is the equivalent generic form.
#
# Plots are PNG, not PDF (per direct user request). The original's `_na_count_plots.pdf` held
# three separate plot() calls as three PDF pages within one file — PNG has no equivalent concept of
# multiple pages in one file, so this one splits into three separate PNGs (na_count_vs_se,
# na_count_vs_score, score_vs_se) rather than being merged into a single combined image.
```

## code/02_data_processing/3_Merging_filtering/bin/build_exclude_list_round1.R

Original header:

```
#!/usr/bin/env Rscript
# build_exclude_list_round1.R — derive the first-round sample exclusion list (low depth + known
# duplicates) from frozen, already-decided inputs.
# Source: 02_data_processing/3_Merging_filtering/01_gencove_axiom_merging_filtering.sh
#         (the R block between QC4 and DA_MERGED_GENCOVE_AXIOM_QC4B)
#
# This is a mechanical derivation, not a human judgment call itself — the judgment already
# happened upstream (the 0.3x depth threshold, and which samples are "duplicates to keep" in the
# frozen file below). Regenerated here (rather than just consumed) as a faithfulness check: the
# result can be diffed against the already-committed data/02_data_processing/list_exclude.txt.
#
# The `fam` join key is `IID`, which — despite this file's `IID` column looking like a short
# "dogID" rather than a raw sequencing sample ID — is verified correct: joining by `IID` exactly
# as written reproduces the original's own documented "454 NA -- axiom mainly" comment precisely.
```

## code/02_data_processing/3_Merging_filtering/bin/build_factor_covariates.R

Original header:

```
#!/usr/bin/env Rscript
# build_factor_covariates.R — per-factor age qcovar + sex covar files
# Source: 02_data_processing/3_Merging_filtering/03_create_files_per_factor_mlma.R
#   F1: lines 93-114, F2: lines 166-187, F3: lines 239-260
#
# `left_join(famQC6, data4, by="IID")` in the original references `data4`, which is never
# defined anywhere in this script or the ones before it in this pipeline (only data1/data5/data6
# exist). `data1` (BUILD_DATA1's output) is the only object in scope with both an `age` and a
# `sex_numeric` column — exactly what this block goes on to use (`famQC6age$age`,
# `famQC6age$sex_numeric`) — and it joins cleanly 1:1 on IID (data1's covariates source has
# exactly one row per dog). Fixed to `data1` here; flagged for confirmation rather than silently
# assumed, since `data4` is undefined and this is the only candidate that makes the rest of the
# block's column references resolve.
#
# `data1$IID` is cast to character in build_data1.R before being written out, but that cast
# doesn't automatically survive a plain-text round trip: read back here with `read.table()`, R
# infers the column's type from its content alone. In practice this pipeline's dogIDs are a mix
# of plain numeric IDs and non-numeric ones (Axiom's hyphenated IDs, GENCOVE's SRR accessions —
# confirmed present in the real data6 checkpoint), so real data1.txt almost certainly infers back
# as character on its own, matching `famQC6$IID` (explicitly cast below) with no issue. A
# synthetic, all-numeric-IID test file *did* reproduce a real crash here (`left_join` errors on
# "incompatible types" under a modern dplyr instead of silently coercing like older versions did)
# — an edge case this defensive, explicit re-cast rules out for free rather than leaving it to
# depend on every batch of real dogIDs happening to include a non-numeric one.
#
# The original's `write.table()` calls for both the qcovar and covar files are commented out for
# all three factors (F1/F2/F3 blocks are identical on this point) — as literally run, the original
# never wrote age*.qcovar/sex*.covar to disk. Written here uncommented per direct user confirmation
# (2026-09-01): these files are the only way to get age/sex covariates to 04_gwas (not yet
# converted), and the process's stated purpose depends on them existing. Several diagnostic prints
# around these blocks (`head(famQC6age.qcovar)`, `table(...$age)`, `head(famQC6sex)`,
# `head/tail(famQC6sex.covar)`, `table(...$sex)`) are also not reproduced in the log below — left
# as-is per the same per-factor-log-asymmetry precedent already accepted in 1_Survey_data; none of
# them affect the written qcovar/covar values.
```

## code/02_data_processing/3_Merging_filtering/bin/build_factor_fam.R

Original header:

```
#!/usr/bin/env Rscript
# build_factor_fam.R — per-factor NA-filtered phenotype fam file
# Source: 02_data_processing/3_Merging_filtering/03_create_files_per_factor_mlma.R
#   F1: lines 45-81, F2: lines 118-156, F3: lines 191-229
#
# Genericized across F1/F2/F3 via arguments. Real per-factor differences preserved:
#   - which item columns feed the NA count (F1: item7/153/154/155, F2: item93/95/150,
#     F3: item145/146/147/148/149/151/152)
#   - the NA-count exclusion threshold (F1: >=3, F2: >=2, F3: >=7) — literal per-factor values
#     from the original, not derived from item count.
#
# `mutate(F1modi=ifelse(na_count>=3,NA,F1))` in the original references a bare `na_count`
# that is never defined anywhere in this script (only `na_count_F1`/`na_count_F2`/`na_count_F3`
# are, one line above each use) — F2 and F3 already use their own `na_count_F{2,3}` correctly;
# F1's is the odd one out. Fixed here to use the per-factor na_count column consistently, per
# explicit confirmation this is a bug to fix (not a literal-reproduction case), since a bare
# `na_count` would be an undefined-object error in a batch run and there is no other candidate
# object it could sensibly refer to.
#
# The fam-file write here was commented out in the original (`#write.table(...)`) even though the
# very next step (a `plink --fam` call) reads this exact file back in — the same
# "write-commented-but-immediately-read-back" idiom already resolved elsewhere in this pipeline
# (see build_data1.R). Written for real here since nothing else in the repo freezes this file.
#
# The original's per-factor NA-count diagnostic (`plot(data6$na_count_F1, ...)`) is a separate,
# uncaptured interactive plot() call that feeds nothing downstream. Captured here as a real PNG
# per factor, per direct user request — one file per factor (not combined across F1/F2/F3).
```

## code/02_data_processing/3_Merging_filtering/bin/build_item_fam.R

Original header:

```
#!/usr/bin/env Rscript
# build_item_fam.R — per-item age-NA-filtered phenotype fam file
# Source: 02_data_processing/3_Merging_filtering/04_create_files_per_item_polmm.sh, lines 17-38
#
# The original shows only one worked example (item155) with the comment "OBS repeat for each
# item, ie. 14 items!!!!!!!!!!!!!!!!!!!!!!" — genericized here across all 14 items via the
# `item_id` argument. Unlike the per-factor fam files, there is no NA-count/multi-item threshold
# step here: each item is filtered only on age-NA, since these are single-item ordinal scores,
# not combined factor scores. This fam-file write was NOT commented out in the original (unlike
# the per-factor equivalent), so it was always a real output.
#
# Validated against the real, frozen data6_2024-10-14.txt checkpoint: item155's raw NA count
# comes out to 728 here, not the 629 the original script's own comment for this exact line
# documents (and 1006 after the age-NA filter, not the comment's 913) — while Stage E's F1/F2/F3
# NA counts (which also depend on data6's item columns, including item155) match their own
# original comments exactly against this same file. This suggests data6_2024-10-14.txt's item155
# column differs slightly from whatever data the archived comment in
# 04_create_files_per_item_polmm.sh was computed against, even though both scripts name the same
# frozen file. Not a bug in this script — the transformation is byte-for-byte the original's
# logic; the mismatch is between a stale inline comment and the file actually committed to the
# repo. Trusting the committed data6 file over the stale comment.
```

## code/02_data_processing/3_Merging_filtering/bin/build_item_grab.R

Original header:

```
#!/usr/bin/env Rscript
# build_item_grab.R — per-item combined phenotype/covariate "grab" file for POLMM
# Source: 02_data_processing/3_Merging_filtering/04_create_files_per_item_polmm.sh, lines 54-72
#
# `sex.x` in the original (`mutate(sex_binary=ifelse(sex.x==1,1,0))`) is dplyr's automatic
# collision suffix: both famQC6 (from the plink fam file) and data6 have a `sex` column, so the
# join produces `sex.x` (famQC6's own plink-coded sex) and `sex.y` (data6's sex) — reproduced
# literally here as `sex.x`, matching the original's own comment "works because no sex NA!!!".
```

## code/02_data_processing/3_Merging_filtering/bin/build_size_fam.R

Original header:

```
#!/usr/bin/env Rscript
# build_size_fam.R — SIZE control-phenotype fam file (run via mlma alongside F1/F2/F3)
# Source: data/04_gwas/getSIZE_pheno_genofiles.r, lines 1-48 (SIZE section only — the file's
# second section, STUCK, lines 54-86, is a separate Darwin's Ark item not reproduced here, per
# direct user instruction, 2026-09-01)
#
# The original re-reads data1 from a round trip through disk
# (`data1 <- read.csv("q121_age_sex_DA_MERGED_GENCOVE_AXIOM_QC5.fam", ...)`, itself written by an
# earlier, also-commented-out `write.table(data1, ...)` two lines above) — the exact same
# famQC5+covariates join `BUILD_DATA1` already performs in this pipeline (confirmed: both use the
# identical covariates file, DarwinsDogs_Q121_height_age_sex_batch_GENCOVE_AXIOM_QC4.tsv). Same
# "write-commented-then-read-back" idiom already resolved elsewhere in this pipeline (see
# build_data1.R) — uses this pipeline's own already-built data1.txt directly instead of a fresh,
# redundant famQC5+covariates join.
#
# Q121A.pheno (Darwin's Ark's own SIZE survey-derived phenotype, data/02_data_processing/
# Q121A.pheno) is exposed as the required `size_pheno_file` param.
#
# Validated against real data: real Q121A.pheno + real covariates + a space-normalized copy of
# the dogID-relabel crosswalk fam (3330 dogs, standing in for a real QC5 fam — no real genotype
# data is in this repo) produces 999 -9-recoded (missing) dogs, leaving 2331 — matching the real,
# already-committed data/04_gwas/SIZE_DA_MERGED_GENCOVE_AXIOM_QC6.fam's 2330 rows almost exactly
# (the real QC5 population is 3316, not 3330, fully accounting for the 1-dog difference).
```

## code/02_data_processing/3_Merging_filtering/bin/depth_histograms.R

Original header:

```
#!/usr/bin/env Rscript
# depth_histograms.R — depth distribution diagnostic plots (all merged dogs vs. the final ALLFAM
# subset)
# Source: 02_data_processing/3_Merging_filtering/02_create_ALLFAM_check_stats.sh (the "create
# stats for the samples used in our analyses" R block, dated later than the rest of the script)
#
# hist_1 uses a `DA_merged` object that isn't defined anywhere in the original script itself — it
# only exists because this was run as a continuation of the same interactive R session as
# 01_gencove_axiom_merging_filtering.sh, where DA_merged (fam_qc3modi joined with the depth file)
# was already in scope. Recomputed fresh here from the same two frozen inputs already used in
# Stage C (BUILD_EXCLUDE_LIST_ROUND1) — a deterministic join on frozen files, so recomputing it
# gives an identical object to reusing it, without adding a new cross-stage output just for one
# diagnostic plot.
#
# `DA_all$IID <- DA_all$sampleID` here has no `as.character()` cast, unlike Stage C's equivalent
# line — verified this doesn't cause a type-mismatch join error (unlike an analogous uncast join
# tried elsewhere): fam_2584dogs.txt (which ALLFAM is restricted to) contains only clean,
# all-numeric dogIDs, so both sides of this join infer as the same type. Not a bug here.
```

## code/02_data_processing/3_Merging_filtering/bin/filter_data6_to_qc5.R

Original header:

```
#!/usr/bin/env Rscript
# filter_data6_to_qc5.R — restrict + reorder data6 to exactly QC5's dogID set, in QC5's own
# row order
#
# BUILD_FACTOR_FAM/BUILD_ITEM_FAM build their per-factor/per-item fam files directly from data6,
# one row per data6 row, which FILTER_FACTOR_QC6/FILTER_ITEM_QC6 then substitute into QC5's
# .bed/.bim via plink's strictly *positional* --fam substitution (same mechanic already
# documented in APPLY_EXCLUDE_LOWPASS12/FILTER_DOGID_RELABEL_FAM). data6 is a frozen checkpoint,
# untouched by and unrelated to this pipeline's own genotype-side sample exclusions (LowPass-12,
# the dog-3094 duplicate pair) — as committed, it still has all 3328 rows of the *original*
# (pre-fix) QC5 population, while QC5 itself is now 3316 rows after those exclusions. Left
# unreconciled, every downstream --fam substitution would fail: the .bed's byte size is fixed by
# QC5's 3316 individuals, and a 3328-row --fam file can't be loaded against it.
#
# Fixed here, per direct user instruction (2026-09-01), by restricting data6 to exactly QC5's own
# fam's dogID list, in QC5's own row order (via match(), not a plain ID filter) — order matters
# as much as membership for a positional substitution; data6's own row order has no reason to
# match QC5's. Errors loudly if any QC5 dog is missing from data6 (would otherwise silently
# produce a malformed downstream fam) rather than filtering leniently and hoping the sets already
# line up.
```

## code/02_data_processing/4_Plotting/bin/plot_pca.R

Original header:

```
#!/usr/bin/env Rscript
# plot_pca.R — PCA scatter plots, colored/shaped by breed
# Source: 02_data_processing/4_Plotting/plot_PCA.R
#
# The original reads breed metadata from a literal placeholder path
# (`/path/to/DarwinsArk_20220715_dogs.csv`) — a real path scrubbed before this script was
# committed. The real file has since been added to the repo
# (data/02_data_processing/DarwinsArk_20220715_dogs.csv, columns confirmed: dog, purebred,
# breed1_inputted) — exposed here as a required `dog_breeds_csv` param rather than hardcoded to
# that path, consistent with how every other frozen data/02_data_processing input is handled
# throughout this project.
#
# The axis labels' variance percentages ("PCA1 (Variance = 37.3%)") are hardcoded in the
# original, baked into the visible plot from whatever historical eigenval file the author had on
# hand. Computed here instead from the real eigenval file, per direct user instruction
# (2026-09-01): a plot's own axis label should never state a variance percentage that doesn't
# match the data actually plotted.
#
# The original's first plot (uncolored, pch only — `col=col.vector` is commented out in its own
# `plot()` call) is dropped entirely here, per direct user instruction (2026-09-01): its legend
# listed breed colors ("other breed"/"Labrador Retriever") that never appeared on the actual
# uncolored plot — a real, pre-existing inconsistency in the original, not worth reproducing.
# Only the second, fully colored/shaped-by-breed plot is kept.
```

## code/02_data_processing/4_Plotting/bin/plot_factor_distr.R

Original header:

```
#!/usr/bin/env Rscript
# plot_factor_distr.R — F1/F2/F3 factor-score distribution histograms, stacked vertically
# Source: 02_data_processing/4_Plotting/plot_factor_distr.R
#
# Genericized across F1/F2/F3 via a small per-factor table — the three original blocks are
# identical except: the title/x-axis label text, the N= annotation, and whether axis.title.x is
# actively styled (only F3's block has it uncommented — the bottom plot in the vertical stack is
# the only one that needs its x-axis title styled). Every other margin/padding value is identical
# across all three and reproduced exactly (see the per-plot user request for exact-margin
# fidelity, 2026-09-01).
#
# The N= annotation and the shared xlim/ylim were hardcoded in the original (baked into the
# visible plot from whatever historical run the author had). Computed here instead from the real
# data, per direct user instruction (2026-09-01) — N via nrow() (each QC6 fam already excludes
# missing-phenotype dogs via plink's --prune, so nrow() is exactly the plotted N), and xlim/ylim
# via the real combined range/max bin height across all three factors (with a small pad), since
# the original applies one shared range across all three plots for visual comparability.
```

## code/02_data_processing/4_Plotting/bin/plot_item_factor_hist.R

Original header:

```
#!/usr/bin/env Rscript
# plot_item_factor_hist.R — factor-score histograms (row) + per-item response histograms (3 rows)
# Source: 02_data_processing/4_Plotting/plot_item&factor_hist.R
#
# Genericized via per-plot config tables (one row per plot) rather than 17 near-duplicated ggplot
# blocks — every functional/styling difference the original actually has (title, x-axis label,
# color, x-axis breaks/limits, title/axis-title margins, y-axis hjust, whether the y-axis is
# blanked) is preserved exactly per plot, per direct user request for exact-margin fidelity
# (2026-09-01), including one genuine oddity reproduced as-is rather than "fixed": item152 has no
# `scale_x_continuous(breaks=1:6)` call unlike its 6 CCDF3 siblings.
#
# item147's original `limits = c(0.5, 6, 5)` (a 3-element vector where ggplot expects 2 — verified
# harmless on its own, ggplot silently uses only the first two elements via ggplot_build()) is
# fixed here to `c(0.5, 6.5)`, per direct user observation: this codebase's other Swedish-locale
# tells (e.g. `SD2VÄRDE` in plot_factor_distr.R) make "6,5" typed as a European decimal comma —
# silently split by R into two vector elements inside c(...) instead of one 6.5 — a far more
# likely explanation than a deliberately-open-ended upper bound.
#
# N= annotations and the shared axis limits (x+y for the 3 factor-score plots — continuous mirt
# scores, whose real range won't match the original's hardcoded historical -1.6..3.5 in general;
# y only for the 14 item plots, since their x-axis is the fixed 1-6 Likert instrument scale, kept
# literal) are computed from the real data here rather than the original's hardcoded historical
# values, per direct user instruction (2026-09-01).
```

## code/02_data_processing/3_Merging_filtering/nextflow.config

```
// Not part of the original script — added per direct user instruction. 01_mapping's own
// README documents these 12 samples as "excluded in the end" in this pipeline (no Darwin's
// Ark dogID, no phenotypes), but 01_gencove_axiom_merging_filtering.sh never actually removes
// them. See RELABEL_AND_EXCLUDE_ROUND3's header comment for how they're removed here instead.
```

## code/02_data_processing/3_Merging_filtering/nextflow_schema.json (`exclude_list_lowpass12` description)

```
12 LowPass samples with no Darwin's Ark dogID and no phenotypes (exclude_12_LowPass_samples.txt, data/01_mapping). Not removed by the original script itself; added per direct user instruction — see RELABEL_AND_EXCLUDE_ROUND3.
```
