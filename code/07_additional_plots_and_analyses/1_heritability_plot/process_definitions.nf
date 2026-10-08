// ============================================================================
// process_definitions.nf — 07. Additional plots and analyses — 1. Heritability plot
//
// Source: plot_FIG_1B_heritabilities.R, moved here; setwd()/hardcoded xlsx path and output
// filename replaced with CLI args. Restyled 2026-10-06 to the published Fig 1b layout (human
// reference lines, labeled items); plotted values unchanged.
// ============================================================================

process PLOT_HERITABILITY {
    label 'process_single'
    container 'community.wave.seqera.io/library/gwas_supplement_plots:e626727dcb23673b'

    input:
    path(supp_tables_xlsx, arity: '1')

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    plot_FIG_1B_heritabilities.R ${supp_tables_xlsx} question_heritabilities_by_factor_boxplot.pdf
    """

    stub:
    """
    touch question_heritabilities_by_factor_boxplot.pdf
    """

    output:
    path("question_heritabilities_by_factor_boxplot.pdf"), emit: pdf
}
