// ============================================================================
// main.nf — 07. Additional plots and analyses — 1. Heritability plot
//
// Fig 1b: per-question heritability boxplot by CCD factor, from the published
// Supplementary Tables (ST1/ST4) — no live Google Sheet dependency.
//
// Source: plot_FIG_1B_heritabilities.R (moved here, see process_definitions.nf header)
// ============================================================================

include { validateParameters } from 'plugin/nf-schema'

include { PLOT_HERITABILITY } from './process_definitions.nf'

workflow {

    main:
    validateParameters()

    ch_supp_tables = channel.fromPath(params.supp_tables_xlsx, checkIfExists: true)

    plot_out = PLOT_HERITABILITY(ch_supp_tables)

    publish:
    pdf = plot_out.pdf

}

// Output publishing mode is symlink by default. This saves space over copying.
output {
    pdf {
        path 'heritability_plot'
    }
}
