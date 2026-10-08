// ============================================================================
// main.nf — 07. Additional plots and analyses — 5. Cross-species transcriptomic heatmap
//
// Extended Data Fig. 3: heatmap of pairwise Spearman rank correlation coefficients between
// human and dog pseudo-bulk expression profiles for the 165 OCD/dogCD network genes, across
// matched cell types.
//
// Source: ../dogCD_ED-fig4_human_dog_heatmap.R (see process_definitions.nf header — collaborator
// script, left in place; bin/cross_species_transcriptomic_heatmap.R is the wired copy). Verified
// 2026-09-30: all 400 cells of the recomputed correlation matrix match the published figure
// exactly — see this directory's README.
// ============================================================================

include { validateParameters } from 'plugin/nf-schema'

include { CROSS_SPECIES_TRANSCRIPTOMIC_HEATMAP } from './process_definitions.nf'

workflow {

    main:
    validateParameters()

    ch_human_dog_rds = channel.fromPath(params.human_dog_rds, checkIfExists: true)
    ch_seed_gene_list = channel.fromPath(params.seed_gene_list, checkIfExists: true)

    heatmap_out = CROSS_SPECIES_TRANSCRIPTOMIC_HEATMAP(ch_human_dog_rds, ch_seed_gene_list)

    publish:
    correlation_by_celltype = heatmap_out.correlation_by_celltype
    correlation_matrix      = heatmap_out.correlation_matrix
    heatmap_plot            = heatmap_out.heatmap_plot

}

// Output publishing mode is symlink by default. This saves space over copying.
output {
    correlation_by_celltype {
        path 'cross_species_transcriptomic_heatmap'
    }
    correlation_matrix {
        path 'cross_species_transcriptomic_heatmap'
    }
    heatmap_plot {
        path 'cross_species_transcriptomic_heatmap'
    }
}
