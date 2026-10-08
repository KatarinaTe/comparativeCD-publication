// ============================================================================
// main.nf — 05. Network — 5. GO semantic clustering workflow
//
// GO semantic clustering (Fig 2c/2e) -> OCD "NA" detail (Fig 2d), fold-enrichment
// barchart (Fig 3b), and the community gene-enrichment tileplots (Fig 3d/e).
//
// Source: bin/{01_semantic_clustering.R, 02_plot_GO_clusters.R, 03_OCD_NA_detail.R,
//          plot_fold_enrichment.R, plot_gene_enrichment_grids.R,
//          plot_gene_enrichment_grids.v2.R, plot_genes_for_multiple_terms.R}
// See process_definitions.nf and archive/conversion_notes/05_network.md
// for the staging/flagging decisions behind each process.
// ============================================================================

include { validateParameters } from 'plugin/nf-schema'

include { SEMANTIC_CLUSTERING            } from './process_definitions.nf'
include { PLOT_GO_CLUSTERS               } from './process_definitions.nf'
include { OCD_NA_DETAIL                  } from './process_definitions.nf'
include { PLOT_SEMANTIC_SIMILARITY       } from './process_definitions.nf'
include { PLOT_FOLD_ENRICHMENT           } from './process_definitions.nf'
include { PLOT_GENE_ENRICHMENT_GRIDS     } from './process_definitions.nf'
include { PLOT_GENE_ENRICHMENT_GRIDS_V2  } from './process_definitions.nf'
include { PLOT_GENES_FOR_MULTIPLE_TERMS  } from './process_definitions.nf'

workflow {

    main:
    validateParameters()

    // The 6 GO enrichment TSVs both 01_semantic_clustering.R and plot_fold_enrichment.R read
    // (via the "../Figure4/..." staging fix — see process_definitions.nf) — .collect() so each
    // single-file channel can be reused across both consumer processes.
    ch_dep_only_go  = channel.fromPath(params.dep_only_go, checkIfExists: true).collect()
    ch_dep_ccd_go   = channel.fromPath(params.dep_ccd_go, checkIfExists: true).collect()
    ch_ocd_only_go  = channel.fromPath(params.ocd_only_go, checkIfExists: true).collect()
    ch_ocd_ccd_go   = channel.fromPath(params.ocd_ccd_go, checkIfExists: true).collect()
    ch_sch_only_go  = channel.fromPath(params.sch_only_go, checkIfExists: true).collect()
    ch_sch_ccd_go   = channel.fromPath(params.sch_ccd_go, checkIfExists: true).collect()
    // dogCD: single-species, dog-only network -- one enrichment file, no run1/run2 pair.
    // Only consumed by SEMANTIC_CLUSTERING (plot_fold_enrichment.R's own script is unchanged
    // and doesn't know about it).
    ch_dogcd_go     = channel.fromPath(params.dogcd_go, checkIfExists: true).collect()

    ch_cache_gpt = channel.fromPath(params.cache_gpt_labels, checkIfExists: true)
    ch_goa_gaf   = channel.fromPath(params.goa_human_gaf, checkIfExists: true).collect()
    ch_seed_csv  = channel.fromPath(params.seed_genes_csv, checkIfExists: true).collect()

    semantic_out = SEMANTIC_CLUSTERING(
        ch_dep_only_go, ch_dep_ccd_go, ch_ocd_only_go, ch_ocd_ccd_go, ch_sch_only_go, ch_sch_ccd_go,
        ch_dogcd_go, ch_cache_gpt, params.openai_api_key
    )

    // Both downstream plotting scripts consume the same per-term output file.
    ch_per_term = semantic_out.per_term.collect()

    go_clusters_out = PLOT_GO_CLUSTERS(ch_per_term)
    ocd_na_out      = OCD_NA_DETAIL(ch_per_term)
    semantic_sim_out = PLOT_SEMANTIC_SIMILARITY(ch_per_term)

    fold_out = PLOT_FOLD_ENRICHMENT(
        ch_dep_only_go, ch_dep_ccd_go, ch_ocd_only_go, ch_ocd_ccd_go, ch_sch_only_go, ch_sch_ccd_go
    )

    // PLOT_GENE_ENRICHMENT_GRIDS, _V2, and PLOT_GENES_FOR_MULTIPLE_TERMS all load this same
    // cache rather than building it themselves — PLOT_FOLD_ENRICHMENT is the actual builder
    // (see process_definitions.nf) — .collect() to reuse it across all 3.
    ch_all_go_rds = fold_out.all_go_rds.collect()

    grids_out = PLOT_GENE_ENRICHMENT_GRIDS(ch_all_go_rds, ch_goa_gaf, ch_seed_csv)

    grids_v2_out = PLOT_GENE_ENRICHMENT_GRIDS_V2(ch_all_go_rds, ch_goa_gaf, ch_seed_csv)
    multi_out    = PLOT_GENES_FOR_MULTIPLE_TERMS(ch_all_go_rds, ch_goa_gaf, ch_seed_csv)

    publish:
    clusters_labeled       = semantic_out.clusters_labeled
    per_term               = semantic_out.per_term
    cutoff_scan             = semantic_out.cutoff_scan
    go_clusters_pdf         = go_clusters_out.pdf
    go_clusters_png         = go_clusters_out.png
    ocd_na_pdf              = ocd_na_out.pdf
    ocd_na_png              = ocd_na_out.png
    semantic_similarity_pdf = semantic_sim_out.pdf
    semantic_similarity_png = semantic_sim_out.png
    semantic_similarity_csv = semantic_sim_out.csv
    fold_enrichment_pdf     = fold_out.pdf
    fold_enrichment_png     = fold_out.png
    fold_enrichment_all_go_rds = fold_out.all_go_rds
    grids_pdf               = grids_out.pdf
    grids_csv               = grids_out.csv
    grids_legend            = grids_out.legend
    grids_overlap           = grids_out.overlap
    grids_v2_pdfs           = grids_v2_out.tileplot_pdfs
    grids_v2_csvs           = grids_v2_out.tileplot_csvs
    grids_v2_overlap        = grids_v2_out.overlap
    multi_pdfs              = multi_out.tileplot_pdfs
    multi_csvs              = multi_out.tileplot_csvs
    multi_overlap           = multi_out.overlap

}

// Output publishing mode is symlink by default. This saves space over copying.
output {
    clusters_labeled {
        path 'semantic_clustering'
    }
    per_term {
        path 'semantic_clustering'
    }
    cutoff_scan {
        path 'semantic_clustering'
    }
    go_clusters_pdf {
        path 'go_clusters_plot'
    }
    go_clusters_png {
        path 'go_clusters_plot'
    }
    ocd_na_pdf {
        path 'ocd_na_detail'
    }
    ocd_na_png {
        path 'ocd_na_detail'
    }
    semantic_similarity_pdf {
        path 'semantic_similarity'
    }
    semantic_similarity_png {
        path 'semantic_similarity'
    }
    semantic_similarity_csv {
        path 'semantic_similarity'
    }
    fold_enrichment_pdf {
        path 'fold_enrichment'
    }
    fold_enrichment_png {
        path 'fold_enrichment'
    }
    fold_enrichment_all_go_rds {
        path 'fold_enrichment'
    }
    grids_pdf {
        path 'gene_enrichment_grids'
    }
    grids_csv {
        path 'gene_enrichment_grids'
    }
    grids_legend {
        path 'gene_enrichment_grids'
    }
    grids_overlap {
        path 'gene_enrichment_grids'
    }
    grids_v2_pdfs {
        path 'tileplots_pdfs_and_csvs'
    }
    grids_v2_csvs {
        path 'tileplots_pdfs_and_csvs'
    }
    grids_v2_overlap {
        path 'tileplots_pdfs_and_csvs'
    }
    multi_pdfs {
        path 'genes_for_multiple_terms'
    }
    multi_csvs {
        path 'genes_for_multiple_terms'
    }
    multi_overlap {
        path 'genes_for_multiple_terms'
    }
}
