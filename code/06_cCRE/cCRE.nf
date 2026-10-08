// ============================================================================
// main.nf — 06. cCRE — emission-state-plot workflow
//
// Chromatin state emission-probability / genome-annotation-enrichment plot.
//
// Source: 06_cCRE/emission_state_UU.R
// ============================================================================

include { validateParameters } from 'plugin/nf-schema'

include { EMISSION_STATE_PLOT } from './process_definitions.nf'

workflow {

    main:
    validateParameters()

    ch_probabilities = channel.fromPath(params.probabilities, checkIfExists: true)
    ch_annotations   = channel.fromPath(params.annotations, checkIfExists: true)

    plot_out = EMISSION_STATE_PLOT(ch_probabilities, ch_annotations)

    publish:
    pdf = plot_out.pdf
    png = plot_out.png

}

// Output publishing mode is symlink by default. This saves space over copying.
output {
    pdf {
        path 'cCRE'
    }
    png {
        path 'cCRE'
    }
}
