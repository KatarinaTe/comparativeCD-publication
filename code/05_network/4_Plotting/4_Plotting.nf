// ============================================================================
// main.nf — 05. Network — 4. Plotting workflow (Fig. 2a-b)
//
// Source: 05_network/4_Plotting/Fig2_network-overlap.R
// ============================================================================

include { validateParameters } from 'plugin/nf-schema'

include { PLOT_VENN_HIERARCHY        } from './process_definitions.nf'
include { ADD_VENN_REGION_LABELS     } from './process_definitions.nf'
include { BUILD_NETWORK_GENE_COUNTS  } from './process_definitions.nf'
include { PLOT_FIG2A                 } from './process_definitions.nf'

workflow {

    main:
    validateParameters()

    ch_ocd_ccd_genes = channel.fromPath(params.ocd_ccd_genes, checkIfExists: true)
    // .collect(): reused below by both PLOT_VENN_HIERARCHY and ADD_VENN_REGION_LABELS -- a queue
    // channel can only be consumed by one process, so this needs to be a value channel to feed two.
    ch_gene_level_network_basis_table = channel.fromPath(params.gene_level_network_basis_table, checkIfExists: true).collect()

    plot_b_out = PLOT_VENN_HIERARCHY(ch_gene_level_network_basis_table)
    venn_regions_out = ADD_VENN_REGION_LABELS(ch_gene_level_network_basis_table)

    counts_out = BUILD_NETWORK_GENE_COUNTS(
        channel.fromPath(params.dogcd_z, checkIfExists: true),
        channel.fromPath(params.dogcd_seeds, checkIfExists: true),
        channel.fromPath(params.dogcd_hier, checkIfExists: true),
        channel.fromPath(params.ocd_z, checkIfExists: true),
        channel.fromPath(params.ocd_seeds, checkIfExists: true),
        channel.fromPath(params.ocd_hier, checkIfExists: true),
        channel.fromPath(params.dep_z, checkIfExists: true),
        channel.fromPath(params.dep_seeds, checkIfExists: true),
        channel.fromPath(params.dep_hier, checkIfExists: true),
        channel.fromPath(params.sch_z, checkIfExists: true),
        channel.fromPath(params.sch_seeds, checkIfExists: true),
        channel.fromPath(params.sch_hier, checkIfExists: true),
        channel.fromPath(params.ocd_ccd_zcomb, checkIfExists: true),
        channel.fromPath(params.ocd_ccd_genes, checkIfExists: true),
        channel.fromPath(params.dep_ccd_zcomb, checkIfExists: true),
        channel.fromPath(params.dep_ccd_hier, checkIfExists: true),
        channel.fromPath(params.sch_ccd_zcomb, checkIfExists: true),
        channel.fromPath(params.sch_ccd_hier, checkIfExists: true),
    )

    plot_a_out = PLOT_FIG2A(counts_out.counts)

    publish:
    figure2b_dogcd_pdf     = plot_b_out.dogcd_pdf
    figure2b_dogcd_png     = plot_b_out.dogcd_png
    figure2b_ocd_pdf       = plot_b_out.ocd_pdf
    figure2b_ocd_png       = plot_b_out.ocd_png
    figure2b_ocd_dogcd_pdf = plot_b_out.ocd_dogcd_pdf
    figure2b_ocd_dogcd_png = plot_b_out.ocd_dogcd_png
    gene_table_with_venn_regions = venn_regions_out.table
    network_gene_counts    = counts_out.counts
    figure2a               = plot_a_out.plot

}

// Output publishing mode is symlink by default. This saves space over copying.
output {
    figure2b_dogcd_pdf {
        path 'fig2'
    }
    figure2b_dogcd_png {
        path 'fig2'
    }
    figure2b_ocd_pdf {
        path 'fig2'
    }
    figure2b_ocd_png {
        path 'fig2'
    }
    figure2b_ocd_dogcd_pdf {
        path 'fig2'
    }
    figure2b_ocd_dogcd_png {
        path 'fig2'
    }
    network_gene_counts {
        path 'fig2'
    }
    figure2a {
        path 'fig2'
    }
    gene_table_with_venn_regions {
        path 'fig2'
    }
}
