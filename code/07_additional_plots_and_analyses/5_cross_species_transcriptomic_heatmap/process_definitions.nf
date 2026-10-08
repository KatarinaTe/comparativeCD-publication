// ============================================================================
// process_definitions.nf — 07. Additional plots and analyses — 5. Cross-species transcriptomic
// heatmap
//
// Source: ../dogCD_ED-fig4_human_dog_heatmap.R (collaborator script, left untouched at its
// original location). bin/cross_species_transcriptomic_heatmap.R is a literal reproduction —
// see its own header comment for exactly what changed and why (readRDS()/gene-list paths as CLI
// args, abbrev_map/cell_abbrevs/col_fun_cor filled in, cell_types narrowed to the 20 plotted
// types, explicit saveRDS()/pdf() calls).
// ============================================================================

process CROSS_SPECIES_TRANSCRIPTOMIC_HEATMAP {
    label 'process_high'
    container 'community.wave.seqera.io/library/cross_species_expression:3b296196c0c71539'

    input:
    path(human_dog_rds, arity: '1')
    path(seed_gene_list, arity: '1')

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    cross_species_transcriptomic_heatmap.R ${human_dog_rds} ${seed_gene_list} .
    """

    stub:
    """
    touch human_dog_correlation_by_celltype.csv
    touch human_dog_correlation_matrix.rds
    touch human_dog_heatmap.pdf
    """

    output:
    path("human_dog_correlation_by_celltype.csv"), emit: correlation_by_celltype
    path("human_dog_correlation_matrix.rds"),      emit: correlation_matrix
    path("human_dog_heatmap.pdf"),                 emit: heatmap_plot
}
