// ============================================================================
// main.nf — 07. Additional plots and analyses — 3. Cross-species expression
//
// Fig 4a-b: dogCD/OCD network gene-set enrichment against the Siletti et al. 2023 human brain
// snRNA-seq atlas, by cell type and by brain region.
//
// Source: network_gene_expression_enrichment_plots_code_updated.r (see process_definitions.nf
// header — collaborator script, left in place; bin/cross_species_expression.R is the wired copy)
// ============================================================================

include { validateParameters } from 'plugin/nf-schema'

include { CROSS_SPECIES_EXPRESSION } from './process_definitions.nf'

workflow {

    main:
    validateParameters()

    ch_siletti_rds = channel.fromPath(params.siletti_rds, checkIfExists: true)

    expr_out = CROSS_SPECIES_EXPRESSION(ch_siletti_rds)

    publish:
    cell_specificity    = expr_out.cell_specificity
    tissue_specificity  = expr_out.tissue_specificity
    celltype_plot       = expr_out.celltype_plot
    region_plot         = expr_out.region_plot

}

// Output publishing mode is symlink by default. This saves space over copying.
output {
    cell_specificity {
        path 'cross_species_expression'
    }
    tissue_specificity {
        path 'cross_species_expression'
    }
    celltype_plot {
        path 'cross_species_expression'
    }
    region_plot {
        path 'cross_species_expression'
    }
}
