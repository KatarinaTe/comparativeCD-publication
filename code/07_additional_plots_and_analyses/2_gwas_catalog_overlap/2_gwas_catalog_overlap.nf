// ============================================================================
// main.nf — 07. Additional plots and analyses — 2. GWAS-catalog overlap
//
// The "12 of 27 dogCD GWAS regions overlap a GWAS-catalog psychiatric gene" claim's mechanical
// half: per-gene psychiatric/OCD/temperament hit-count join, reproduced from local repo-tracked
// inputs only (no live Google Sheets access needed).
//
// Source: reproduce_gwas_catalog_overlap.R (moved here, see process_definitions.nf header)
// ============================================================================

include { validateParameters } from 'plugin/nf-schema'

include { REPRODUCE_GWAS_CATALOG_OVERLAP } from './process_definitions.nf'

workflow {

    main:
    validateParameters()

    ch_catalog_hits     = channel.fromPath(params.catalog_hits_file, checkIfExists: true)
    ch_name_conversion  = channel.fromPath(params.name_conversion_file, checkIfExists: true)
    ch_frozen_overlap   = channel.fromPath(params.frozen_overlap_file, checkIfExists: true)
    ch_frozen_gene_list = channel.fromPath(params.frozen_gene_list_file, checkIfExists: true)

    overlap_out = REPRODUCE_GWAS_CATALOG_OVERLAP(
        ch_catalog_hits, ch_name_conversion, ch_frozen_overlap, ch_frozen_gene_list
    )

    publish:
    signif_hits          = overlap_out.signif_hits
    hit_summary           = overlap_out.hit_summary
    gene_list_with_hits    = overlap_out.gene_list_with_hits
    crosscheck             = overlap_out.crosscheck

}

// Output publishing mode is symlink by default. This saves space over copying.
output {
    signif_hits {
        path 'gwas_catalog_overlap'
    }
    hit_summary {
        path 'gwas_catalog_overlap'
    }
    gene_list_with_hits {
        path 'gwas_catalog_overlap'
    }
    crosscheck {
        path 'gwas_catalog_overlap'
    }
}
