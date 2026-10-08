// ============================================================================
// process_definitions.nf — 05. Network — cCRE
//
// Source: 06_cCRE/emission_state_UU.R
//
// New container: r-ggplot2/r-patchwork/r-dplyr/r-tidyr/r-scales, all conda-forge, no exotic deps.
// See env/06_cCRE/README.md for the build story and why the R implementation (not the
// duplicate Python notebook, emission_state_UU.ipynb) was chosen — ggsave's width=16/height=7/
// dpi=300 exactly matches the deposited emission_states_plot.png's pixel dimensions (4800x2100),
// confirming R is what actually produced it; the Python duplicate renders at different dimensions.
// ============================================================================

process EMISSION_STATE_PLOT {
    label 'process_single'
    container 'community.wave.seqera.io/library/emission_state:0737d6956e5a4a7d'

    input:
    path(probabilities, arity: '1')
    path(annotations, arity: '1')

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    emission_state_plot.R ${probabilities} ${annotations} emission_states_plot.pdf emission_states_plot.png
    """

    stub:
    """
    touch emission_states_plot.pdf
    touch emission_states_plot.png
    """

    output:
    path("emission_states_plot.pdf"), emit: pdf
    path("emission_states_plot.png"), emit: png
}
