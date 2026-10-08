# Archive

This folder holds material that is no longer part of the active analysis but is kept for
provenance rather than deleted outright.

- `conversion_notes/` — per-stage notes on the decisions made while converting each stage's
  original scripts into the containerized Nextflow pipelines under `code/`. This content used to
  live inline in stage `README.md` files and script headers; it has been moved here so those files
  only document the science (what a pipeline does, its inputs/outputs, how to run it) rather than
  the process of building it.
- `legacy_readmes/` — pre-conversion `README` files that described the original, non-Nextflow
  scripts for a stage and are now superseded by that stage's own pipeline `README.md`.
- `elinors_analysis_superseded/` — files from a collaborator's original working directory that are
  draft/scratch/duplicate/orphaned relative to what actually feeds a figure, table, or claim in the
  manuscript (see `data/Manuscript_draft/REVIEW_AND_SUGGESTIONS.md` for the specific figure-by-figure
  mapping and reasoning behind each move). Content is unchanged from the original — only its
  location moved. (Named before this repo's later naming-by-content convention below; not yet
  renamed — see the note in `go_semantic_clustering_originals/` for why.)
- `go_semantic_clustering_originals/` — the rest of that same working directory's GO-semantic-
  clustering/fold-enrichment/community-tileplot material: pre-generated figure outputs (the exact
  PDFs/PNGs/CSVs behind Fig 2c-e/3b/3d-e, superseded once Nextflow regenerates them) and misc
  project/review files. The 9 files actually needed to reproduce the paper's results (7 scripts, 2
  frozen inputs) moved instead to `code/05_network/5_GO_semantic_clustering/bin/` and
  `data/05_network/5_GO_semantic_clustering/` — those are the only ones still outside archive, per
  direct instruction that only Nextflow-related files needed to produce the paper's results should
  live outside `archive/`. (`elinors_analysis_superseded/` above predates this folder's split-out and
  rename and hasn't been revisited yet.)
- `gwas_catalog_enrichment_originals/` — the same working directory's separate GWAS-catalog "12 of
  27" claim: its exploratory scripts and the raw per-gene `gwas_cache/` query cache. Kept apart from
  `go_semantic_clustering_originals/` since it's an unrelated analysis, not a superseded draft of it.
  The 3 files actually needed to reproduce the paper's results (`gwas_catalog_hits.csv`,
  `name_conversion_v2.tsv`, `gene_list_gwas_catalog_overlap_with_region.csv`) moved to
  `data/07_additional_plots_and_analyses/` (2026-09-21) — same "pull the real inputs out, archive
  the rest" treatment `go_semantic_clustering_originals/` already got.
- `cCRE_fig4_debugging_aids/` (2026-09-22) — two one-off diagnostic scripts written during the
  2026-09-17 Fig. 4d/4e debugging session: `diagnose_other_tissues.sh` (per-tissue chrom/bp
  diagnostic used to isolate the EpicDog mammary-gland pooling bug) and
  `compare_UU_epicdog_elements_bedtools.sh` (a pure-bedtools reimplementation of Fig. 4d's R/valr
  comparison, written to check whether a numeric discrepancy was a valr-vs-bedtools tooling issue —
  it wasn't). Cluster cross-checks used to find/verify bugs, not part of the final pipeline;
  confirmed no other file references either script before moving them.
- `testing_CRE_overlap_superseded/` — the original, non-pipeline scripts behind Fig. 4e/4f
  (`testing_CRE_overlap.R`, its successor `testing_CRE_overlap_260617.R`, and the two bedtools
  recipes `check_CRE_overlap.sh`/`SIZE_check_UU_CRE_overlap.sh`). `testing_CRE_overlap.R` was
  archived 2026-09-18 as a confirmed-superseded duplicate of `_260617.R`. The other three joined
  it 2026-09-21 once `_260617.R`'s real content (the logic actually behind the two published
  panels, most of the rest being exploratory/superseded duplicate work) was converted into
  `code/06_cCRE/cre_overlap/` — see that stage's README for exactly what was kept vs. dropped.
- `network_gene_expression_enrichment_superseded/` — the pre-2026-09-15 version of
  `code/07_additional_plots_and_analyses/network_gene_expression_enrichment_plots_code_updated.r`
  (Fig 4a-b, Siletti et al. 2023 human snRNA-seq cross-species comparison), superseded by the
  collaborator's own updated version: it fixes a real bug (the neuron-only per-brain-region GSEA
  used `human_tissue_specificity_scores`, computed from a `human_neurons` subset that was never
  actually created) and adds the Figshare provenance comment for the one deposited input,
  `Siletti_human_subset_normalised.rds`.
- `go_term_comparison_unused/` (2026-09-21) — the full, working Nextflow-converted
  `GO_TERM_COMPARISON` workflow (`code/`) and its run config (`analyses/`), formerly
  `code/05_network/4_Plotting` and `analyses/05_network/4_Plotting`. Unlike everything else in
  this folder, this isn't superseded or a draft — it's real, tested, runnable code. Moved here
  per KatarinaTe's direct instruction because its output (GO-term "run1 single-species vs. run2
  cross-species" difference volcano plots, for OCD/DEP/SCH) isn't tied to any figure/table in the
  current manuscript, and she intends to remove it from the repo entirely. Confirmed before moving
  that nothing else depends on it (see `REPRODUCIBILITY_AUDIT.md`'s 2026-09-21 update).
  `Fig2_network-overlap.R`, the other, unrelated script that used to share the same directory
  (Fig. 2a-b, genuinely used), stayed in place at `code/05_network/4_Plotting/`.
- `fig2b_euler_diagram_superseded/` (2026-10-01) — `plot_network_overlap.R`, the Fig. 2b
  stacked-bar + Euler diagram process, superseded by area-proportional Venn diagrams
  (`PLOT_VENN_HIERARCHY` / `bin/plot_venn_hierarchy.R`) after KatarinaTe provided an updated Fig. 2.
  Kept for provenance, including the real "wrong gene list" bug found and fixed in it 2026-09-22
  (see its own header comment and `code/05_network/4_Plotting/README.md`'s Notes).
- `go_semantic_clustering_pre_dogcd/` (2026-10-03) — the pre-dogCD version of
  `01_semantic_clustering.R` and its matching `cache_gpt_labels.top_community.cut60.txt` (140
  clusters), superseded when the collaborator's revised script added dogCD as a 4th network
  (needed to feed the corrected Fig. 2e semantic-similarity panel) plus per-cluster cohesion
  metrics and a cutoff-sensitivity scan. See `code/05_network/5_GO_semantic_clustering/README.md`'s
  Notes for the full story, including the one new cluster (`cluster_141`) this produced.
- `fig2e_similarity_draft_superseded/` (2026-10-03) — the earlier Fig. 2e semantic-similarity
  reconstructions (`fig2e_similarity_heatmap_DRAFT*`, `*_with_dogCD*`, `*_v3_full_legend*`, the
  draft markdown doc, and this folder's own README) built from the manuscript Results text alone,
  before the collaborator's real figure legend and intermediate values became available. Superseded
  by `code/05_network/5_GO_semantic_clustering`'s `PLOT_SEMANTIC_SIMILARITY` once that was verified
  against the collaborator's own computed values. The bug that made every draft attempt here wrong
  on the OCD–Depression cell: GOSemSim's `combine="BMA"` is size-weighted (pools all terms from both
  sets before averaging), not the equal-weighted average of the two directional means the real
  method uses — only visible for badly size-imbalanced pairs like OCD (8 terms) vs. Depression
  (129). Confirmed by manually reproducing the collaborator's own intermediate numbers exactly. See
  `code/05_network/5_GO_semantic_clustering/README.md`'s Notes for the full writeup.
- `supplementary_tables_superseded/` (2026-10-06) — earlier drafts of the manuscript's
  Supplementary Tables workbook, superseded by `data/Manuscript_draft/NATURE_SupplementaryTables_261006.xlsx`
  (the current version, which `1_heritability_plot` reads). `NATURE_SupplementaryTables.xlsx` is
  the previously tracked copy (last edited 2026-09-29 09:07), `NATURE_SupplementaryTables_260929.xlsx`
  a later same-day draft, and `NATURE_SupplementaryTables_261005.xlsx` the 2026-10-05 draft
  (originally saved as `NATURE_SupplementaryTables (1).xlsx`; renamed by date on archiving). Content
  unchanged. ST4 heritabilities are identical across all four versions; only the 2026-10-05 and
  later versions carry ST1's `Short name` column, which Fig 1b's item labels need.
- `cCRE_overlap_tables_superseded/` (2026-10-06) — three bp-overlap tables that used to sit in
  `data/06_cCRE/` (`epicdog_tissue_bp_overlaps.txt`, `brain_region_bp_overlaps.txt`,
  `SIZE_brain_region_bp_overlaps.txt`), built 2026-09-17/21 from pre-final GWAS regions. No pipeline
  reads them; `06c-cre-overlap` regenerates all three from the final regions, and its live output
  (2026-09-23 run) reproduces Supplementary Table 13 and the Results text's per-region odds ratios
  exactly, whereas these older copies don't (e.g. EpicDog dogCD OR 1.08 vs. the published 1.06).
  Content unchanged.
