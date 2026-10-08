// ============================================================================
// process_definitions.nf — 07. Additional plots and analyses — 3. Cross-species expression
//
// Source: network_gene_expression_enrichment_plots_code_updated.r (collaborator script, left
// untouched at its original location). bin/cross_species_expression.R is a literal
// reproduction — only the readRDS() path replaced with a CLI arg, and explicit
// saveRDS()/ggsave() calls added for the 2 derived matrices and 2 plots the original built
// in-memory but never wrote to disk.
// ============================================================================

process CROSS_SPECIES_EXPRESSION {
    label 'process_high'
    container 'community.wave.seqera.io/library/cross_species_expression:3b296196c0c71539'

    input:
    path(siletti_rds, arity: '1')

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    cross_species_expression.R ${siletti_rds}
    """

    stub:
    """
    touch human_cell_specificity_scores.rds
    touch human_tissue_specificity_scores.rds
    touch human_network_gsea_results_plot.pdf
    touch human_network_per_region_gsea_results_plot.pdf
    """

    output:
    path("human_cell_specificity_scores.rds"),               emit: cell_specificity
    path("human_tissue_specificity_scores.rds"),              emit: tissue_specificity
    path("human_network_gsea_results_plot.pdf"),               emit: celltype_plot
    path("human_network_per_region_gsea_results_plot.pdf"),   emit: region_plot
}
