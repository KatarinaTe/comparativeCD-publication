# Conversion notes — 04_gwas

Decisions and discrepancies noted while converting this stage's original scripts into the
Nextflow pipelines under `code/04_gwas/`. Moved out of the `.nf`/`process_definitions.nf` header
comments and `bin/*` script docstrings so those files document only the science/pipeline itself.

## code/04_gwas/1_mlma-loco/1_mlma-loco.nf

**GRM/REML removed again (2026-09-15), reversing the 2026-09-11 reinstatement below.**
KatarinaTe's PR comment gave the actual citation chain for Fig 1b: GitHub repo
(`VistaSohrab/dog-gwas-heritability-nextflow`, calculation) → Supplementary Table 4 (published
result) → `plot_FIG_1B_heritabilities.R` (plotting) — confirming this repo's own GRM/REML
computation was never meant to feed Fig 1b and is redundant with it. `BUILD_GRM`/`RUN_REML`/
`PARSE_HERITABILITY` and `bin/parse_heritability.R` are gone again;
`plot_FIG_1B_heritabilities.R` was repointed at `NATURE_SupplementaryTables.xlsx`'s `ST4`/`ST1`
sheets directly (confirmed to have the exact columns/layout the script already expected) instead
of the private Google Sheet it read from before. See the entries below for the reinstatement this
reverses.

```
Per-phenotype mlma-loco GWAS, for the 3 CCD factors and SIZE (the control phenotype) — then
plink clumping + sumstats export for all 4 phenotypes (see process_definitions.nf's header for
why SIZE is included in clumping/export despite no literal source for that step).

SIZE is kept in its own published subdirectories from gene-finding (clumping) onward — the
clumping step is where genes first get attached to a signal (--clump-range), and the sumstats
files feed downstream finemapping/region-storage — per direct user instruction (2026-09-02) not
to mix SIZE's regions in with the real factors'. Raw mlma-loco GWAS output (before gene-finding
starts) stays in one shared directory across all 4 phenotypes.

**GRM/REML reinstated (2026-09-11), reversing the 2026-09-01 removal below.** The manuscript's
Fig 1b (SNP-based heritability boxplot,
code/07_additional_plots_and_analyses/plot_FIG_1B_heritabilities.R) needs these exact numbers and
currently gets them from a private, undocumented Google Sheet instead. Per direct user
instruction, GRM-building + both REML calls are back, for all 4 phenotypes (F1/F2/F3 + SIZE),
feeding a new heritability_estimates.csv with columns matching what
plot_FIG_1B_heritabilities.R already expects (phenotype, sample_size,
GREML_LDS_constrained_heritability, GREML_LDS_no_constraint_heritability), so that script could
eventually be pointed at this file instead (not done here — out of scope, plot_FIG_1B_heritabilities.R
itself is untouched).

**SIZE included in REML — flagged, not obviously right.** This stage's established pattern (see
above) is SIZE runs through everything the factors do, so BUILD_GRM/RUN_REML/PARSE_HERITABILITY
run for SIZE too. But unlike mlma-loco/clumping/sumstats-export, there's no established-by-instruction
precedent for this specific case: MLMA_PER_FACTOR.sh's own REML lines only ever mention "Factor
1... factor 2-3", never SIZE, and plot_FIG_1B_heritabilities.R itself explicitly filters its
question-level data to phenotype_name %in% c("CCDF1","CCDF2","CCDF3") — it never reads a SIZE row
at all. Included anyway for pipeline uniformity and because it's cheap and harmless (an extra CSV
row nothing downstream currently reads); flagging in case a SIZE heritability estimate is
specifically unwanted.

**"GREML_LDS_*" column names vs. plain GREML commands — a preexisting naming mismatch, not
introduced here.** MLMA_PER_FACTOR.sh's own REML `--out` values are literally named
`*.REML.no-lds` (no LD-score adjustment) — no `--reml-lds` flag or any LD-score step appears
anywhere in this script. Yet plot_FIG_1B_heritabilities.R's Google Sheet columns are named
"GREML LDS constrained heritability" / "GREML LDS noConstraint heritability". This pipeline's new
heritability_estimates.csv reproduces the *values* MLMA_PER_FACTOR.sh's commands actually produce
(plain GREML, not LD-score-stratified GREML) under the "LDS"-named columns purely to match
plot_FIG_1B_heritabilities.R's expected shape for eventual drop-in use — it does not resolve
whether "LDS" in the manuscript's column names is a misnomer, or refers to a different LD-score
step run elsewhere and never captured in this script. Flagged, not decided.

**BUILD_GRM omits --thread-num, unlike RUN_MLMA_LOCO/RUN_REML.** MLMA_PER_FACTOR.sh's own
`--make-grm` call (line 16) never passes --thread-num at all (only the REML calls do, explicitly,
at 16) — reproduced literally, so BUILD_GRM defaults to gcta64's own thread default (likely 1)
rather than task.cpus. This could make genome-wide GRM-building slower than necessary; flagged as
a literal-fidelity choice rather than silently adding a flag the original never had.
```

## code/04_gwas/1_mlma-loco/process_definitions.nf

```
MLMA_PER_FACTOR.sh's own example is F1 only ("do the same for factor 2-3"), genericized here
across F1/F2/F3 via a config channel — plus SIZE, run through the identical mlma-loco step per
direct user confirmation (2026-09-01: SIZE's real historical run used the same command shape
and covariates as the factors).

Clumping and sumstats export ALSO run for SIZE, per direct user instruction (2026-09-02): SIZE
should come out the other end as clumped regions (after downstream finemapping) just like the
factors, even though clump_plink_factors.sh/export_gwas_sumstats_CCD_MLMA.R only ever name
F1/F2/F3 — same category of confirmed extension-by-instruction as SIZE's own mlma-loco run
above, not a literal source. CLUMP_FACTOR/EXPORT_SUMSTATS branch their output prefix on
pheno_id ('SIZE' vs. "CCD${pheno_id}") rather than hardcoding "CCD", matching
1_mlma-loco.nf's own bfilePrefix() convention for the same SIZE-has-no-"CCD"-prefix quirk.

The original's GRM-building step (line 16) and both REML/heritability calls (lines 19-39) were
deliberately not reproduced initially, per direct user instruction (2026-09-01) — go straight from
the QC6 genotypes to mlma-loco, on the reasoning that --mlma-loco computes its own per-chromosome
(leave-one-chromosome-out) GRMs internally and never takes a precomputed GRM as input (the
original's own commented-out line 45, `#gcta64 --mlma --grm ...`, is the non-LOCO variant that
*would* use it — not what's actually run). **Reinstated 2026-09-11** (see the entry above in this
file): that reasoning about --mlma-loco not consuming a precomputed GRM still holds — BUILD_GRM's
GRM is consumed only by RUN_REML, never by RUN_MLMA_LOCO — but the REML/heritability numbers
themselves turned out not to be dead work: the manuscript's Fig 1b needs them, and was getting
them from an undocumented Google Sheet instead. BUILD_GRM/RUN_REML/PARSE_HERITABILITY reproduce
lines 16-39 for all 4 phenotypes (F1/F2/F3 + SIZE, see the SIZE-inclusion flag above), one process
each: BUILD_GRM (`--make-grm`), RUN_REML (both `--reml` calls back-to-back, sharing one GRM
input), PARSE_HERITABILITY (a new bin/parse_heritability.R parsing both `.hsq` outputs per
phenotype into one combined heritability_estimates.csv row, via Nextflow's
`collectFile(keepHeader: true)` rather than an extra combining process).

RUN_MLMA_LOCO's --out uses the full QC6 bfile prefix (e.g. CCDF1_DA_MERGED_GENCOVE_AXIOM_QC6_LOCO)
rather than the original's abbreviated CCDF1_QC6_LOCO — cosmetic only (nothing parses this
filename), kept per direct user confirmation (2026-09-01) since it's unambiguous across all 4
phenotypes in a flat published output directory.

`--thread-num` uses `task.cpus` rather than the original's hardcoded 16 — the original's value
was tied to a specific SLURM node allocation (`#SBATCH -n 1`, a full node), not an algorithmic
requirement; thread count doesn't affect GWAS results, only wall-clock time.

Reuses 3_Merging_filtering's container (gcta64 + R/dplyr) for mlma-loco/sumstats export, and
the existing plink container for clumping — no new build needed.
```

## code/04_gwas/2_polmm/process_definitions.nf

```
POLMMgrab_PER_ITEM.sh's own example is item155 only ("OBS!!! Here example for item155: Change
and run for each item"), genericized here across all 14 items via a config channel.

A real inconsistency in the original, noted but not something this conversion could reproduce
either way: POLMMgrab_PER_ITEM.sh's own awk step writes `modi_simuMarkerOutput_...` (lowercase),
but clump_plink_items.sh reads `MODI_simuMarkerOutput_...` (uppercase) — a case mismatch between
the two original scripts. Doesn't arise here since RUN_POLMM_ITEM's real output file is passed
to CLUMP_ITEM directly via a Nextflow channel, not re-derived from a hardcoded filename string.

New container (community.wave.seqera.io/library/polmm_grab) — the GRAB R package (POLMM
method) isn't part of 3_Merging_filtering's existing R stack. Clumping reuses the existing
plink container; sumstats export's R script needs no packages beyond base R.

UNVERIFIED-VERSION CAVEAT: pinned to GRAB 0.2.5, the newest available (conda-forge only goes
back to 0.2.2; CRAN's archive only to 0.2.1). The real analysis may have used an earlier
version (possibly 0.1.1, per the user) that predates both archives and whose GitHub repo
couldn't be located to check. `GRAB.NullModel`/`GRAB.Marker` (the two functions this pipeline
calls) exist under those exact names in 0.2.5, but their defaults/algorithm details in an
actual 0.1.1 are unverified — same category of open caveat as 1_Survey_data's unverified
`mirt()` seeding question. Flagged per direct user confirmation (2026-09-01) rather than
silently assumed equivalent.
```

## code/04_gwas/3_finemap_susie/3_finemap_susie.nf

```
Builds one shared genofile (split per chromosome) that all phenotypes finemap against, then
harmonizes + munges each of the 18 GWAS'd phenotypes' sumstats (3 factors + 14 items, from the
original scripts, plus SIZE — see below), then runs PolyFun's finemapper.py for every one of the
113 windows in assets/finemap_regions.csv (104 hand-picked) + assets/finemap_regions_size.csv
(9, mechanically extracted from SIZE's own already-run results, added 2026-09-02). SIZE goes
through the identical processes but is kept in its own published subdirectories throughout
(sumstats/size, finemap/size) — same "keep SIZE separated" split as 1_mlma-loco's own clumping.

[...]

The 17 phenotypes harmonized/munged/finemapped in the original scripts (STUCK excluded
throughout this stage — see process_definitions.nf's header). SIZE is wired in separately below
— no literal source names it here either, but data/04_gwas/SIZE_260603_CSnonfunct.xlsx (added
2026-09-02) is direct evidence SIZE WAS finemapped in the real analysis, with its own region
windows recoverable from that file (see extract_size_regions.py).
```

## code/04_gwas/3_finemap_susie/process_definitions.nf

```
Structurally different from 1_mlma-loco/2_polmm: instead of one worked example genericized
across a phenotype/item list, the finemapper.py calls are ~100 individual, hand-picked
(phenotype, chr, start, end) windows accumulated over an iterative, exploratory process — some
loci re-run at 2-4 progressively narrower widths. Per direct user instruction (2026-09-02):
reproduce every one of these calls literally (not just the apparently-final/narrowest one per
locus). The full table was extracted MECHANICALLY (a one-off parsing script, not hand-typed)
from the 4 original shell scripts into assets/finemap_regions.csv (104 rows) to avoid
transcription error — see that file's own header.

STUCK is excluded from harmonize/munge/finemap entirely, matching the literal original: no
finemap-region windows exist for it anywhere, and it was excluded from GWAS entirely per
3_Merging_filtering's own scope. Its deposited fam file (STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.fam,
referenced below) is used only as a sample-inclusion mask for the one shared genofile all
phenotypes finemap against — not as an analysed phenotype.

SIZE, by contrast, IS harmonized/munged/finemapped here (added 2026-09-02) — despite
2a_run_harmonize_sumstats_CCD.sbatch's own SIZE line being commented out in the original. No
finemap-region windows for SIZE exist in any of the 4 region-table scripts either (there's no
generic "same shape as the factors" command for finemapping the way there was for clumping — a
finemap window is a human-picked choice per locus). But data/04_gwas/SIZE_260603_CSnonfunct.xlsx
(added to `main` 2026-09-02, after this pipeline's initial implementation) is direct evidence
SIZE WAS finemapped in the real analysis — its own already-run region filenames double as the
missing window list, mechanically extracted into assets/finemap_regions_size.csv (9 rows) by
extract_size_regions.py, same "extract mechanically, don't hand-type" practice as the 104-row
table. See 3_finemap_susie.nf's header and extract_size_regions.py's own header for the full
reasoning (including why the credible-set boundary columns in that spreadsheet are NOT used as
the finemap window — they're the credible set's own narrower bounds, not the region analysed to
find it).

New container required: PolyFun (github.com/omerwe/polyfun) is an external tool with no
versioned releases (rolling git history only) — not vendored anywhere in this repo, not
pip/conda-installable as a package. UNVERIFIED-VERSION CAVEAT: pinned to commit 3e657a1066
(2024-07-28), the latest commit that existed before the original analysis's own recorded run
date (its working directory is literally named ".../2024-10-11/") — a traceable anchor, unlike
GRAB's unverified-version caveat in 2_polmm, but still not a confirmed exact match, since the
analyst could have used an older clone. Built from PolyFun's own polyfun.yml (conda-forge only,
python=3.8, r-susier==0.11.92 via rpy2) plus a `git clone` + `git checkout` layered on top (see
env/04_gwas/3_finemap_susie/polyfun.yml and env/04_gwas/README.md). Reused for harmonize/munge/
finemap alike, since harmonize only needs pandas (already in PolyFun's own env).

The genofile-creation R snippet (1_finemap_create_genofiles.sh) is itself informal/incomplete —
its own header says these are notes on "these scripts use these geno inputs ... created like
these", and it references `famQC5` without ever showing it being read. Reconstructed using
3_Merging_filtering's own published QC5 fam (DA_MERGED_GENCOVE_AXIOM_QC5.fam) — see
bin/build_combined_fam.R's header.

The original's `--hwe 'midp' 0.00000000000000000001` argument order (modifier before the
p-value) is reproduced literally below despite plink 1.9's documented syntax being
`--hwe <p-value> [midp]` (value first) — unverified whether plink actually accepts this reversed
order or silently misparses it; flagged here rather than silently reordered, no real data
available yet to test which way it actually behaves.
```

## code/04_gwas/4_qc_signals_IGV/4_qc_signals_IGV.nf

```
Source: 04_gwas/4_qc_signals_IGV/check_signals_igv_manually.sh (personal scratch notes, not a
finished script — see process_definitions.nf's header for the full rationale).
```

## code/04_gwas/4_qc_signals_IGV/process_definitions.nf

```
Source: 04_gwas/4_qc_signals_IGV/check_signals_igv_manually.sh — personal "NOTE TO SELF"
scratch-work (its own text: "this analysis needs to be redone... these are now out of date"),
not a finished script to reproduce literally. This automates everything leading up to the one
genuinely manual step (a human visually judging IGV read pileups) — which SNP to check per
credible set, which BAM-having individuals to load, and re-generating a fresh, non-overlapping
batch of individuals if an earlier batch had no usable reads at that exact position.

Two real ID/data gaps were found and resolved while building this, without touching
01_mapping/3_Merging_filtering's own code or published outputs:

1) The lead-SNP list doesn't need 5_gwas_regions (on hold) at all for this purpose: PolyFun's
   finemapper.py (SuSiE) already tags every SNP with a CREDIBLE_SET column in 3_finemap_susie's
   own finemap_results output (0 = not in a set; verified by reading finemapper.py's source).

2) 01_mapping's BAM filenames use the raw sequencing sampleID; 3_Merging_filtering's QC5/QC6
   genotypes use "dogID" (relabeled via a positional fam substitution with no explicit lookup
   table published anywhere). DA_MERGED_GENCOVE_AXIOM_QC3modi.fam — already-deposited raw data,
   predates the relabel — turns out to carry BOTH the original sampleID and the final dogID
   (IID) columns in the same row; verified row-for-row against the relabel fam's own IID
   sequence. That file is reused directly as the crosswalk (see bin/build_bam_dogid_keep_list.py's
   header for the full reasoning).
```

## env/04_gwas/README.md

```
**Unverified-version caveat:** pinned to 0.2.5, the newest version available on conda-forge (which
only goes back to 0.2.2) and CRAN (archive only to 0.2.1). The real analysis this pipeline
reproduces may have used an earlier version (possibly 0.1.1) that predates both archives — its
GitHub repo couldn't be located to check defaults/algorithm details against. Flagged per direct
user confirmation (2026-09-01) rather than assumed equivalent; see
`code/04_gwas/2_polmm/process_definitions.nf`'s header comment.

The A1/A2-swap concern specifically (does the export script's relabeling still match 0.2.5's
allele-frequency/beta convention?) **is** confirmed, by real execution: `GRAB.Marker`'s output
labels its `Info` field `CHR:POS:REF:ALT` and reports `AltFreq`/`beta` relative to ALT, matching
the assumption the original script's swap depends on (2026-09-02). This narrows, but doesn't close,
the broader unverified-version caveat above — see `bin/export_gwas_sumstats_polmm.R`'s header.
```

```
own `polyfun.yml` conda spec is layered with a `git clone` + `git checkout` build step, following
this project's own documented convention for source-installed tools
(see `~/.claude/skills/docs/decisions/0015-conda-run-command-needs-run-prefix.md`).
```
(This path reference was deleted outright, not archived for reuse — it points to a local
AI-assistant tooling file with no documentary value to a reader of this repo. Recorded here only
for completeness of what was removed.)

```
**Unverified-version caveat, with a firmer anchor than `2_polmm`'s:** PolyFun has no tags/releases,
so pinned to commit `3e657a1066` (2024-07-28) — the latest commit that existed before the original
analysis's own recorded run date (its own working directory is literally named `.../2024-10-11/`).
This is a traceable anchor, not a guess, but still not a confirmed exact match — the analyst could
have cloned an older commit. Flagged per this project's established unverified-version-caveat
practice (same category as `2_polmm`'s GRAB pin and `1_Survey_data`'s unverified `mirt()` seeding);
see `code/04_gwas/3_finemap_susie/process_definitions.nf`'s header comment.
```

## code/04_gwas/1_mlma-loco/bin/export_gwas_sumstats_mlma.R

```
N is computed here from the real .phen file's row count (the actual number of dogs the GWAS
was run on) rather than the original's hardcoded historical value, per the same
compute-real-values-not-stale-hardcoded-ones instruction already applied in 4_Plotting
(2026-09-01) — this is baked into the exported sumstats file itself, not just a log/comment.
```

## code/04_gwas/1_mlma-loco/bin/parse_heritability.R

```
New script (2026-09-11), added for the GRM/REML reinstatement — see the entries above under
1_mlma-loco.nf and process_definitions.nf for why. Sample size (n) is read from the .hsq file
itself (GCTA reports the n it actually used for that REML fit) rather than recomputed from the
.phen file the way export_gwas_sumstats_mlma.R computes N — simpler, and no reason to expect
GCTA's n to differ from the .phen row count here.
```

## code/04_gwas/1_mlma-loco/README.md

```
Created (2026-09-11) — no per-substage README previously existed for any 04_gwas substage (only
the stage-wide, pre-conversion code/04_gwas/README_04.gwas.md), unlike 02_data_processing's
per-substage READMEs. Written from scratch, in that established per-substage template
(Pipeline Overview / Input Data / Output Data / Process Map / Notes), covering this substage's
whole current pipeline (mlma-loco, GRM/REML heritability, clumping, sumstats export) rather than
only the new REML pieces, since there was no existing file to append a new-outputs section to.
```

## code/04_gwas/2_polmm/bin/export_gwas_sumstats_polmm.R

```
The A1/A2 swap (`A1 = wood$A2, A2 = wood$A1`) is intentional, not a bug — reproduced exactly per
the original's own heavily-flagged comment: "OBS! BETA must be for the A1 allele!!! ... changing
the A2 to become A1". GRAB/POLMM's own A1/A2 convention doesn't match what this sumstats format
(and downstream finemapping) expects BETA to be measured against, so the original author
deliberately relabels the columns here.

N (= AltCounts/AltFreq) is already computed from real per-marker values in the original, not a
hardcoded historical constant like the mlma-loco export's original had — no fix needed here.

Confirmed against a real run of GRAB 0.2.5 (the pinned, unverified-version container): GRAB's own
GRAB.Marker output labels its Info field "CHR:POS:REF:ALT", and reports AltFreq/AltCounts/beta/
seBeta relative to ALT (the second allele) — but reformat_polmm_output.sh's positional split
labels REF as "A1" and ALT as "A2" (matching Info's REF-then-ALT order literally, not by meaning).
So before this swap, freq/BETA are actually measured against A2, not A1 — exactly the mismatch
the original author's comment describes. The swap makes the final A1 the ALT allele, consistent
with freq/BETA both being reported against it. Verified end-to-end with GRAB's own bundled
example data (2026-09-02).
```

## code/04_gwas/2_polmm/bin/reformat_polmm_output.sh

```
Column 1 is just an echo of the input .bim file's marker ID (GRAB doesn't compute it) — the
"chr:pos" assumption only holds because both upstream genotype pipelines
(01_mapping/2_Axiom_imputation) set variant IDs via `bcftools annotate --set-id '%CHROM:%POS'`
before merging. Verified against GRAB's own bundled generic example data (plain "SNP_N" marker
names, not chr:pos) that this awk step would silently produce an empty Position column if that
upstream convention were ever missing — confirmed real risk, but not one this pipeline
introduces or needs to guard against, since the convention is already established upstream.
```

## code/04_gwas/2_polmm/bin/run_polmm_item.R

```
Source: 04_gwas/2_polmm/POLMMgrab_PER_ITEM.sh, lines 12-39 (item155's own worked example,
"OBS!!! Here example for item155: Change and run for each item")

Genericized across all 14 items via the item_id argument. Step 1 fits the null model on the
pruned genotype (EigenInput, ~1M SNPs, matching the original's own comment); step 2 runs the
full-genome marker test using that fitted model.
```
(Only the source-line citation and "matching the original's own comment" clause were trimmed;
the rest of this docstring was already plain present-tense documentation and was kept.)

## code/04_gwas/3_finemap_susie/bin/build_combined_fam.R

```
The original snippet reads QC5's own fam into `famQC5` but never shows that read (its own header
comment says these lines just document "these scripts use these geno inputs ... created like
these" — provenance notes, not a clean standalone script). Reconstructed here as the first
argument, using 3_Merging_filtering's own plain, whole-cohort published QC5 fam
(DA_MERGED_GENCOVE_AXIOM_QC5.fam) rather than one of that stage's several phenotype-specific QC5
fam variants (SIZE_..._QC5.fam, CCD${factor}_..._QC5.fam, CCDitem${id}_..._QC5.fam — same QC5
sample set, each with a different phenotype column substituted in). This is the only one of them
that makes sense here: this script builds ONE shared genofile used by every phenotype's
finemapping, so a phenotype-specific starting fam would be an arbitrary, unmotivated pick among
equally-eligible candidates — and whichever variant was used, only its FID/IID/F/M/sex columns
would matter anyway, since the phenotype column gets overwritten immediately below (all QC5-stage
fam variants describe the same underlying sample set, just with different phenotype values). Not
verified against real data (no real QC5 genotype data available yet, same long-standing blocker
as 2_Axiom_imputation/1_mlma-loco) — flagging only that this specific inference is unverified,
not that the choice itself is ambiguous.
```

## code/04_gwas/4_qc_signals_IGV/bin/build_bam_dogid_keep_list.py

```
DA_MERGED_GENCOVE_AXIOM_QC3modi.fam (a raw, already-deposited data file — QC3 predates the
relabel) turns out to already carry both columns directly: `sampleID` and `IID` (the dogID),
row-for-row consistent with the relabel fam (verified: their IID/dogID sequences match exactly).
So it doubles as the crosswalk this needs, with zero changes to 3_Merging_filtering's own code or
published outputs required.
```

## code/04_gwas/4_qc_signals_IGV/bin/extract_credible_set_snps.py

```
PolyFun's finemapper.py (SuSiE method) already tags every SNP with a `CREDIBLE_SET` column in its
own output (0 = not in any credible set; 1, 2, ... = which credible set it belongs to) — see
finemapper.py's `finemap()` method (SuSiE branch), which builds this directly from the fitted
susieR object's own `sets` attribute. So "which SNPs are worth checking in IGV" doesn't need the
still-on-hold 5_gwas_regions collaborator scripts at all for this narrower purpose — it's already
sitting in every one of 3_finemap_susie's ~104 result files.

[...] `P` survives unrenamed into the finemap output from 3_finemap_susie's own sumstats export
(`munge_polyfun_sumstats.py` normalizes whatever p-value column name it finds to exactly `P` —
confirmed by reading its source — and finemapper.py's SuSiE branch copies the full sumstats row
through into its own output, so `P` is just carried along).
```

## code/04_gwas/4_qc_signals_IGV/bin/select_igv_individuals.py

```
Source: 04_gwas/4_qc_signals_IGV/check_signals_igv_manually.sh's own commented approach ("run for
all 1200 ind for which we have bam-files for the lead SNPs of all associated loci") — but that
script only eyeballs a plink --recode .ped file by hand ("grep A ...ped"); this replaces that
ad-hoc step with an actual genotype-class selection, since the original note is personal, informal
scratch-work (its own header: "NOTE TO SELF"), not a finished script to reproduce literally.

[...]

Input: a plink `--recode A` output (.raw) — NOT the original's `--recode` (.ped), which is
letter-coded and only practical for by-eye grep.
```

(Note: `code/04_gwas/2_polmm/2_polmm.nf` and `code/04_gwas/3_finemap_susie/bin/harmonize_sumstats.py`
contained no conversion-process narration to remove — left unedited.)

## code/04_gwas/5_gwas_regions — built 2026-09-11, removed entirely 2026-09-15

```
dog_gws_regions.sh's own "NOTE TO SELF" placeholder named two scripts as the real analysis that
weren't in the repo yet: cCRE_enrichment_Epic2tissues_Random10times.R and
FishercomparisonForSigRegioPipSNP.R (both per-SNP fine-mapped-cCRE analyses — "Fine-mapped dogCD
SNPs in cCREs" in the manuscript). Both arrived later, committed directly to
code/04_gwas/5_gwas_regions/ (commit 8da62f9), and were built into a 5-process Nextflow pipeline
(EPIC_RANDOM_TISSUE_ENRICHMENT, EPIC_EMPIRICAL_PVALUE, UU_PAIRWISE_ENRICHMENT,
SIGREGIONS_RANDOM_TISSUE_ENRICHMENT, SIGREGIONS_FISHER_COMPARISON), stub-tested, over two frozen
per-SNP chromatin-state-call input tables (AllSnps_StateEnrichment_PromEnhancer_unique.txt,
SigRegionsUniqueSNP_02june2026.txt).

Per direct instruction (2026-09-15): this whole analysis belongs to a different collaborator's own
repository, not this one. Removed entirely — code/04_gwas/5_gwas_regions/, its analyses/ and env/
counterparts, and the two frozen data files — not archived, since it's moving elsewhere rather than
being superseded. dog_gws_regions.sh (the original placeholder) went with it, since its sole purpose
was naming these two scripts.

This does not affect the separate, region-level UU/EPIC dogCD-vs-SIZE cCRE comparisons in
code/06_cCRE/testing_CRE_overlap_260617.R (behind "DogCD GWAS regions are enriched in brain
regulatory elements" / "cCREs in dogCD GWAS regions enriched across brain regions") — a different
script, different methodology (per-region bp-overlap vs. per-SNP tissue calls), predating and
unrelated to this removed pipeline or the PlotsToShare material it arrived alongside (also removed,
2026-09-15 — see the PR comment thread for that decision).

Detail on what was built (process split, the two excluded-tail bugs found in
cCRE_enrichment_Epic2tissues_Random10times.R, the deferred SigRegions empirical-p-value block, the
Cerebellum column-collision data caveat) is not preserved here beyond this summary, since the code
itself no longer exists in this repo to cross-reference against.
```
