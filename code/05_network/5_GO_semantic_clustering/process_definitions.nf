// ============================================================================
// process_definitions.nf — 07. Additional plots and analyses — GO semantic
// clustering / fold-enrichment / community-tileplot chain
//
// Every process below runs one of the collaborator's own, unedited scripts, now living in this
// stage's own bin/ (never modified — read-only collaborator source, same policy as
// every other stage's original-author material; only their location, execute bit, and
// a prepended `#!/usr/bin/env Rscript` shebang have changed, per direct instruction, so
// they can be invoked bare like every other stage's bin/ scripts, resolved via
// Nextflow's automatic bin/-on-PATH — no content/logic was touched).
//
// These scripts were written to be run by hand from a working directory the script
// itself sits in: they use hardcoded relative paths (some of them broken, see below),
// no CLI arguments, and write fixed output filenames. Rather than rewrite their logic,
// each process here stages the exact other inputs each script expects at the exact
// relative path it expects them.
//
// See archive/conversion_notes/05_network.md for the full
// reasoning behind every staging/flagging decision summarised in the short notes below.
// ============================================================================

// ----------------------------------------------------------------------------
// 01_semantic_clustering.R
//
// Updated 2026-10-03 to the collaborator's revised script: adds dogCD (single-species,
// dog-only network — one enrichment file, no run1/run2 pair; its own `human_disease` value,
// `run2` left at 0) as a 4th network alongside DEP/OCD/SCH, needed to feed the corrected Fig. 2e
// (see PLOT_SEMANTIC_SIMILARITY below). Also adds per-cluster cohesion metrics
// (within_cluster_mean_sim/min_sim, representative term) and a clustering-cutoff sensitivity
// scan (silhouette width, cophenetic correlation across cutoffs 0.4-0.8). Until the 2026-10-06
// edit below, content otherwise
// unmodified from the collaborator's own script — same policy as every other stage's
// original-author material; only location, execute bit, and a prepended
// `#!/usr/bin/env Rscript` shebang changed. Verified byte-for-byte reproduction of the
// collaborator's own output files against the real deposited data before integrating — see
// archive/conversion_notes/05_network.md.
//
// Edited 2026-10-06 with KatarinaTe's explicit authorization, for publication: the script's
// hardcoded input directory (a collaborator's Google Drive mount point) now reads the 6
// human-disease GO enrichment TSVs plus the dogCD one from its working directory, where they
// are staged below; and its final two googlesheets4::sheet_write() calls (mirroring the
// per-cluster table and cutoff scan to a shared Google Sheet; they always failed here, with no
// OAuth token in the container) were removed. Both tables are still written to files.
// ----------------------------------------------------------------------------
process SEMANTIC_CLUSTERING {
    label 'process_high'
    container 'community.wave.seqera.io/library/go_semantic_plots:69629d5c597a400b'

    input:
    path(dep_only_go, arity: '1')
    path(dep_ccd_go, arity: '1')
    path(ocd_only_go, arity: '1')
    path(ocd_ccd_go, arity: '1')
    path(sch_only_go, arity: '1')
    path(sch_ccd_go, arity: '1')
    path(dogcd_go, arity: '1')
    path(cache_gpt, arity: '1')
    val(openai_api_key)

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mkdir -p run
    for f in ${dep_only_go} ${dep_ccd_go} ${ocd_only_go} ${ocd_ccd_go} ${sch_only_go} ${sch_ccd_go} ${dogcd_go} ${cache_gpt}; do
        ln -s "../\$f" "run/\$f"
    done

    ( cd run && OPENAI_API_KEY="${openai_api_key}" 01_semantic_clustering.R )
    """

    stub:
    """
    mkdir -p run
    touch "run/GO_clusters_labeled.cut60.txt"
    touch "run/GO_with_without_dogCD.semantic_clustering.top_community.cut60.txt"
    touch "run/GO_cluster_cutoff_scan.top_community.txt"
    """

    output:
    path("run/GO_clusters_labeled.cut60.txt"), emit: clusters_labeled
    path("run/GO_with_without_dogCD.semantic_clustering.top_community.cut60.txt"), emit: per_term
    path("run/GO_cluster_cutoff_scan.top_community.txt"), emit: cutoff_scan
}


// ----------------------------------------------------------------------------
// 02_plot_GO_clusters.R (Fig 2c/2d) — consumes SEMANTIC_CLUSTERING's per-term output
// directly from the working directory (its own list.files(".", ...) glob); no path
// fix needed here. (CORRECTION 2026-10-03: this header previously also credited this
// script with Fig 2e -- it doesn't touch it. 2e is a separate computation, now
// PLOT_SEMANTIC_SIMILARITY below.)
// ----------------------------------------------------------------------------
process PLOT_GO_CLUSTERS {
    label 'process_single'
    container 'community.wave.seqera.io/library/go_semantic_plots:69629d5c597a400b'

    input:
    path(per_term, arity: '1')

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    02_plot_GO_clusters.R
    """

    stub:
    """
    touch "GO_run1_vs_run2_by_semantic_cluster.top_community.cut60.sig5.pdf"
    touch "GO_run1_vs_run2_by_semantic_cluster.top_community.cut60.sig5.png"
    """

    output:
    path("GO_run1_vs_run2_by_semantic_cluster.top_community.cut60.sig5.pdf"), emit: pdf
    path("GO_run1_vs_run2_by_semantic_cluster.top_community.cut60.sig5.png"), emit: png
}


// ----------------------------------------------------------------------------
// 03_OCD_NA_detail.R (Fig 2d) — same per-term input as PLOT_GO_CLUSTERS, same glob.
// ----------------------------------------------------------------------------
process OCD_NA_DETAIL {
    label 'process_single'
    container 'community.wave.seqera.io/library/go_semantic_plots:69629d5c597a400b'

    input:
    path(per_term, arity: '1')

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    03_OCD_NA_detail.R
    """

    stub:
    """
    touch "OCD_NA_detail.top_community.cut60.pdf"
    touch "OCD_NA_detail.top_community.cut60.png"
    """

    output:
    path("OCD_NA_detail.top_community.cut60.pdf"), emit: pdf
    path("OCD_NA_detail.top_community.cut60.png"), emit: png
}


// ----------------------------------------------------------------------------
// plot_semantic_similarity.R (Fig 2e) — added 2026-10-03. NOT a collaborator original: built
// and verified in-session against the collaborator's real figure legend and her own
// intermediate computed values (all 15 pairwise raw-similarity values among DEP/OCD/SCH
// run1/run2 matched her numbers to 3 decimals; the chance-corrected score matched within
// permutation noise) -- see this directory's README.md Notes for the full writeup, including
// the real bug this replaces (GOSemSim's own `combine="BMA"` is size-weighted, not the
// equal-weighted two-direction average the real method uses -- only visible for badly
// size-imbalanced pairs like OCD vs. Depression). Consumes SEMANTIC_CLUSTERING's per-term
// output directly, same as PLOT_GO_CLUSTERS/OCD_NA_DETAIL.
//
// Runtime: 10 pairwise comparisons x 1000 permutations each -- the long pole in this stage
// alongside SEMANTIC_CLUSTERING itself.
// ----------------------------------------------------------------------------
process PLOT_SEMANTIC_SIMILARITY {
    label 'process_high'
    container 'community.wave.seqera.io/library/go_semantic_plots:69629d5c597a400b'

    input:
    path(per_term, arity: '1')

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    plot_semantic_similarity.R ${per_term} fig2e_semantic_similarity
    """

    stub:
    """
    touch fig2e_semantic_similarity.pdf
    touch fig2e_semantic_similarity.png
    touch fig2e_semantic_similarity_results.csv
    """

    output:
    path("fig2e_semantic_similarity.pdf"), emit: pdf
    path("fig2e_semantic_similarity.png"), emit: png
    path("fig2e_semantic_similarity_results.csv"), emit: csv
}


// ----------------------------------------------------------------------------
// plot_fold_enrichment.R (Fig 3b) — reads the same broken "../go_enrichment_tsvs/..." path as
// 01_semantic_clustering.R; same run/ + go_enrichment_tsvs/ staging fix, independently (this
// script never reads 01's cache or outputs).
//
// CORRECTION (2026-09-15, 2nd pass): this is the actual builder of
// all_go_combined.rds — a prior investigation (both the original conversion notes and
// an earlier pass today) wrongly placed that role on plot_gene_enrichment_grids.R,
// which only ever loads it. plot_fold_enrichment.R computes `cache_is_stale` (true
// whenever the file doesn't already exist, which it never does in a fresh task
// directory) and, on that branch, reads the 6 go_enrichment_tsvs/ files and `saveRDS`s the result
// to "all_go_combined.rds" in its own working directory — as a side effect of
// building its own barchart. Emitted here so PLOT_GENE_ENRICHMENT_GRIDS,
// PLOT_GENE_ENRICHMENT_GRIDS_V2, and PLOT_GENES_FOR_MULTIPLE_TERMS (none of which
// build it themselves) can consume it, restoring the dependency chain the original
// conversion notes described but attributed to the wrong process.
// ----------------------------------------------------------------------------
process PLOT_FOLD_ENRICHMENT {
    label 'process_medium'
    container 'community.wave.seqera.io/library/go_semantic_plots:69629d5c597a400b'

    input:
    path(dep_only_go, arity: '1')
    path(dep_ccd_go, arity: '1')
    path(ocd_only_go, arity: '1')
    path(ocd_ccd_go, arity: '1')
    path(sch_only_go, arity: '1')
    path(sch_ccd_go, arity: '1')

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mkdir -p go_enrichment_tsvs run
    ln -s "../${dep_only_go}" "go_enrichment_tsvs/${dep_only_go}"
    ln -s "../${dep_ccd_go}" "go_enrichment_tsvs/${dep_ccd_go}"
    ln -s "../${ocd_only_go}" "go_enrichment_tsvs/${ocd_only_go}"
    ln -s "../${ocd_ccd_go}" "go_enrichment_tsvs/${ocd_ccd_go}"
    ln -s "../${sch_only_go}" "go_enrichment_tsvs/${sch_only_go}"
    ln -s "../${sch_ccd_go}" "go_enrichment_tsvs/${sch_ccd_go}"
    ( cd run && plot_fold_enrichment.R )
    """

    stub:
    """
    mkdir -p run
    touch "run/all_go_combined.rds"
    touch "run/GO_enrichment_OCD_CCD_facets.p5.top10.pdf"
    touch "run/GO_enrichment_OCD_CCD_facets.p5.top10.png"
    """

    output:
    path("run/all_go_combined.rds"), emit: all_go_rds
    path("run/GO_enrichment_OCD_CCD_facets.p5.top10.pdf"), emit: pdf
    path("run/GO_enrichment_OCD_CCD_facets.p5.top10.png"), emit: png
}


// ----------------------------------------------------------------------------
// plot_gene_enrichment_grids.R (Fig 3d/e, "all communities" combined panel) —
// internally still headed "plot_gene_enrichment_grids.v4.R" (a stale docstring, not
// edited). CORRECTION (2026-09-15): this script does NOT build all_go_combined.rds —
// like the other 2 group-4 scripts, it only ever does `stopifnot(file.exists(
// cache_path)); readRDS(cache_path)`. PLOT_FOLD_ENRICHMENT is the actual builder (see
// its own header comment) — its all_go_rds output feeds this process's
// all_go_combined_rds input directly.
//
// Reads seed genes from gene_list_gwas_catalog_overlap_with_region.csv, a frozen export
// of the collaborator's seed-gene Google Sheet tab (switched from a live
// googlesheets4::read_sheet() call 2026-10-06, per KatarinaTe; the frozen file's dog/human
// seed flags and dogCD p-values were checked identical to the live sheet first). Downloads
// goa_human.gaf.gz if absent — staged here from the frozen, already-checked-in copy instead.
//
// Produces no PNG (only a PDF/CSV/legend-info + the overlap TSV), despite a .png
// sitting alongside the other pre-generated outputs in archive/go_semantic_clustering_originals/ —
// flagged, not fabricated (see conversion notes).
// ----------------------------------------------------------------------------
process PLOT_GENE_ENRICHMENT_GRIDS {
    label 'process_high'
    container 'community.wave.seqera.io/library/go_semantic_plots:69629d5c597a400b'

    input:
    path(all_go_combined_rds, arity: '1')
    path(goa_gaf, arity: '1')
    path(seed_genes_csv, arity: '1', stageAs: 'gene_list_gwas_catalog_overlap_with_region.csv')

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mkdir -p run
    ln -s "../${all_go_combined_rds}" "run/all_go_combined.rds"
    ln -s "../${goa_gaf}" "run/${goa_gaf}"
    ln -s "../${seed_genes_csv}" "run/${seed_genes_csv}"
    ( cd run && plot_gene_enrichment_grids.R )
    """

    stub:
    """
    mkdir -p run
    touch "run/OCD_CCD_all_communities_p4_maxts1000_nesting_tileplot.pdf"
    touch "run/OCD_CCD_all_communities_p4_maxts1000_nesting_tileplot_data.csv"
    touch "run/OCD_CCD_all_communities_p4_maxts1000_nesting_tileplot_legend.info.txt"
    touch "run/OCD_CCD_community_gene_overlap.tsv"
    """

    output:
    path("run/OCD_CCD_all_communities_p4_maxts1000_nesting_tileplot.pdf"), emit: pdf
    path("run/OCD_CCD_all_communities_p4_maxts1000_nesting_tileplot_data.csv"), emit: csv
    path("run/OCD_CCD_all_communities_p4_maxts1000_nesting_tileplot_legend.info.txt"), emit: legend
    path("run/OCD_CCD_community_gene_overlap.tsv"), emit: overlap
}


// ----------------------------------------------------------------------------
// plot_gene_enrichment_grids.v2.R (Fig 3d/e, tileplots_pdfs_and_csvs/ per-community
// panels) — internally headed "go_nesting_tileplot.R". Loads
// PLOT_GENE_ENRICHMENT_GRIDS's all_go_combined.rds directly (never builds it), so no
// go_enrichment_tsvs staging is needed here at all.
//
// Its community list (7: C184/C185/C186/C187/C189/C193/C197) and per-community
// output-filename shape ("<ds>_<community>_p<X>_maxts<Y>_nesting_tileplot.*") match
// the files already checked into tileplots_pdfs_and_csvs/, but this script's CURRENT
// pvalue_max (1e-5, -> "p5") does not match those files' "p4" — see conversion notes;
// output filenames are globbed, not hardcoded to "p4", so this process reports
// whatever the script's current config actually computes rather than a guessed value.
// ----------------------------------------------------------------------------
process PLOT_GENE_ENRICHMENT_GRIDS_V2 {
    label 'process_high'
    container 'community.wave.seqera.io/library/go_semantic_plots:69629d5c597a400b'

    input:
    path(all_go_rds, arity: '1')
    path(goa_gaf, arity: '1')
    path(seed_genes_csv, arity: '1', stageAs: 'gene_list_gwas_catalog_overlap_with_region.csv')

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    plot_gene_enrichment_grids.v2.R
    """

    stub:
    """
    touch "OCD_CCD_C185_p5_maxts1000_nesting_tileplot.pdf"
    touch "OCD_CCD_C185_p5_maxts1000_nesting_tileplot_data.csv"
    touch "OCD_CCD_community_gene_overlap.tsv"
    """

    output:
    path("OCD_CCD_*_nesting_tileplot.pdf"), emit: tileplot_pdfs
    path("OCD_CCD_*_nesting_tileplot_data.csv"), emit: tileplot_csvs
    path("OCD_CCD_community_gene_overlap.tsv"), emit: overlap
}


// ----------------------------------------------------------------------------
// plot_genes_for_multiple_terms.R — a 3rd variant of the same "go_nesting_tileplot.R"
// script (own communities: C184/C185/C186/C193/C197; own p_cutoff/max_termsize; no
// "_maxts" in its own output stem). Also loads all_go_combined.rds directly. Named
// explicitly in scope by this stage's task description, but not referenced by any
// figure mapping found in REVIEW_AND_SUGGESTIONS.md or REPRODUCIBILITY_AUDIT.md, and
// its own output-filename shape does not match any already-checked-in deliverable —
// see conversion notes for the full comparison against .v2.R.
// ----------------------------------------------------------------------------
process PLOT_GENES_FOR_MULTIPLE_TERMS {
    label 'process_high'
    container 'community.wave.seqera.io/library/go_semantic_plots:69629d5c597a400b'

    input:
    path(all_go_rds, arity: '1')
    path(goa_gaf, arity: '1')
    path(seed_genes_csv, arity: '1', stageAs: 'gene_list_gwas_catalog_overlap_with_region.csv')

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    plot_genes_for_multiple_terms.R
    """

    stub:
    """
    touch "OCD_CCD_C185_p5_nesting_tileplot.pdf"
    touch "OCD_CCD_C185_p5_nesting_tileplot_data.csv"
    touch "OCD_CCD_community_gene_overlap.tsv"
    """

    output:
    path("OCD_CCD_*_nesting_tileplot.pdf"), emit: tileplot_pdfs
    path("OCD_CCD_*_nesting_tileplot_data.csv"), emit: tileplot_csvs
    path("OCD_CCD_community_gene_overlap.tsv"), emit: overlap
}
