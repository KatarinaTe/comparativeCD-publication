# Conversion notes — 07_additional_plots_and_analyses

This directory now only holds standalone scripts that aren't part of any Nextflow stage yet:
`plot_FIG_1B_heritabilities.R`, `plot_gwas_results.R`, `reproduce_gwas_catalog_overlap.R`, and
`network_gene_expression_enrichment_plots_code.r`. (The GO-semantic-clustering/fold-enrichment/
tileplot Nextflow pipeline that used to live here was moved to
`code/05_network/5_GO_semantic_clustering` on 2026-09-15 — see `archive/conversion_notes/
05_network.md` for its conversion notes.)

## `plot_gwas_results.R` — generalized per direct instruction (2026-09-11 PR comment)

Unlike the GO-semantic-clustering scripts, this file is not read-only original-author content
left in place — KatarinaTe explicitly asked (PR #4 review comment) to "adapt to fit all input
phenotypes and file names," since the version she added only covered SIZE (mlma) and item153
(polmm). Changes made directly to this file rather than left alongside an original:

- Removed the duplicated Manhattan/QQ-plotting code (two near-identical ~100-line blocks) into one
  `plot_gwas_results()` function called once per phenotype/item.
- Fixed a bare `R` on its own line (a leftover interactive-console marker between the polmm and mlma
  blocks) that would error if the original file were run as a script.
- Removed the hardcoded absolute `setwd("/proj/snic2022-6-229/...")` — matches every other script in
  this repo's convention of taking the working directory as given rather than hardcoding a path.
- MLMA-LOCO side now loops over `CCDF1, CCDF2, CCDF3, SIZE`, matching file names to
  `1_mlma-loco/process_definitions.nf`'s real output pattern
  (`${qc6_bed.baseName}_LOCO.loco.mlma`) via `Sys.glob()`, and warns/skips rather than erroring if a
  phenotype's file isn't present (per KatarinaTe: "we don't right now provide SIZE plots anywhere but
  I guess it can stay").
- POLMM side kept parameterized over `polmm_item_ids` (a plain vector to edit) rather than hardcoding
  item153, since POLMM runs per survey item in `2_polmm/process_definitions.nf`, not per CCD factor —
  there's no fixed "all POLMM items" list specified anywhere to enumerate automatically.
- Not wired into a Nextflow process (same as `plot_FIG_1B_heritabilities.R` elsewhere in this
  directory) — no container/stub test exists for it yet.
