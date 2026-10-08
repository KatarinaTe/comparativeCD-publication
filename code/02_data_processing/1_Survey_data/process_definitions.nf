// ============================================================================
// process_definitions.nf — 1. Survey data (EFA/IRT)
//
// Source: 02_data_processing/1_Survey_data/Darwin_dog-CD_3-factorsolution.R
//
// Every diagnostic print and plot in the original script is reproduced as a real output file
// here — sink()-captured logs and png()-captured plots. Plots are PNG rather than the original's
// PDF; RUN_IRT_FACTOR's na_count_plots (originally 3 pages in one PDF) splits into 3 separate
// PNGs accordingly.
// ============================================================================

process RUN_EFA {
    // EFA diagnostics (nfactors/parallel analysis/nScree, KMO, Bartlett, determinant, cronbach's
    // alpha) and fit. None of this feeds the F1/F2/F3 IRT models below — the item groupings they
    // use are a fixed, pre-decided grouping, not derived programmatically from this fit.
    tag   "cohort"
    label 'process_medium'
    container 'community.wave.seqera.io/library/survey_data_efa_irt:f4e0f4b9755aa129'

    input:
    path(response_df, arity: '1')
    path(questions_csv, arity: '1')

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    run_efa.R ${response_df} ${questions_csv}
    """

    stub:
    """
    touch efa_diagnostics.log
    touch parallel_analysis_plot.png
    touch nscree_plot.png
    touch uniqueness_plot.png
    """

    output:
    path("efa_diagnostics.log"),         emit: log
    path("parallel_analysis_plot.png"),  emit: parallel_analysis_plot
    path("nscree_plot.png"),             emit: nscree_plot
    path("uniqueness_plot.png"),         emit: uniqueness_plot
    tuple val("${task.process}"), val('r-base'), eval("Rscript --version 2>&1 | sed 's/.*version //;s/ .*//'"), emit: versions, topic: versions
}


process RUN_IRT_FACTOR {
    // Per-factor IRT fit and factor scores, genericized via arguments rather than three
    // near-duplicate processes. See bin/run_irt_factor.R's header for the real per-factor
    // differences preserved (M2() `type`, which fscores() method feeds the saved CSV, F1's one
    // extra diagnostic plot, and the tracePlot/testInfoPlot/itemInfoPlot theta ranges + color
    // palette).
    tag   "${factor_id}"
    label 'process_medium'
    container 'community.wave.seqera.io/library/survey_data_efa_irt:f4e0f4b9755aa129'

    input:
    tuple val(factor_id), val(item_cols), val(item_weights), val(m2_type), val(score_method),
          val(trace_theta_range), val(testinfo_theta_range), val(iteminfo_theta_range), val(iteminfo_palette)
    path(response_df, arity: '1')

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    run_irt_factor.R ${response_df} ${factor_id} '${item_cols}' '${item_weights}' '${m2_type}' '${score_method}' \\
        '${trace_theta_range}' '${testinfo_theta_range}' '${iteminfo_theta_range}' '${iteminfo_palette}'
    """

    stub:
    """
    touch ${factor_id}_CCD3F.txt
    touch ${factor_id}_diagnostics.log
    touch ${factor_id}_trace_plot.png
    touch ${factor_id}_iteminfo_facet_plot.png
    touch ${factor_id}_testinfo_plot.png
    touch ${factor_id}_iteminfo_plot.png
    touch ${factor_id}_eap_map_plot.png
    touch ${factor_id}_na_count_vs_se_plot.png
    touch ${factor_id}_na_count_vs_score_plot.png
    touch ${factor_id}_score_vs_se_plot.png
    """

    output:
    tuple val(factor_id), path("${factor_id}_CCD3F.txt"),           emit: scores
    path("${factor_id}_diagnostics.log"),                           emit: log
    path("${factor_id}_trace_plot.png"),                            emit: trace_plot
    path("${factor_id}_iteminfo_facet_plot.png"),                   emit: iteminfo_facet_plot
    path("${factor_id}_testinfo_plot.png"),                         emit: testinfo_plot
    path("${factor_id}_iteminfo_plot.png"),                         emit: iteminfo_plot
    path("${factor_id}_eap_map_plot.png"),                          emit: eap_map_plot
    path("${factor_id}_na_count_vs_se_plot.png"),                   emit: na_count_vs_se_plot
    path("${factor_id}_na_count_vs_score_plot.png"),                emit: na_count_vs_score_plot
    path("${factor_id}_score_vs_se_plot.png"),                      emit: score_vs_se_plot
    tuple val("${task.process}"), val('mirt'), eval("Rscript -e \"cat(as.character(packageVersion('mirt')))\""), emit: versions, topic: versions
}


process PLOT_EXTRA_ITEM_SE {
    // F1's one extra diagnostic scatter plot (an item's response vs. factor SE), that F2/F3
    // don't have. Split out from RUN_IRT_FACTOR and invoked only for F1 (see the main workflow)
    // rather than forcing RUN_IRT_FACTOR to declare a conditionally-produced output — Nextflow
    // output paths require a minimum arity of 1, so "produced only for some factors" can't be
    // expressed as an optional output on a process all three factors share. Reads the factor's
    // already-saved CCD3F.txt rather than duplicating any computation.
    tag   "${factor_id}"
    label 'process_single'
    container 'community.wave.seqera.io/library/survey_data_efa_irt:f4e0f4b9755aa129'

    input:
    tuple val(factor_id), path(scores_csv), val(item_col)

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    plot_extra_item_se.R ${scores_csv} ${factor_id} ${item_col}
    """

    stub:
    """
    touch ${factor_id}_extra_item_se_plot.png
    """

    output:
    path("${factor_id}_extra_item_se_plot.png"), emit: plot
}


process PLOT_FACTOR_DISTRIBUTIONS {
    // Final per-factor score distribution plots. Needs all three factors' CCD3F files together —
    // a barrier after the per-factor fan-out.
    tag   "cohort"
    label 'process_single'
    container 'community.wave.seqera.io/library/survey_data_efa_irt:f4e0f4b9755aa129'

    input:
    path(scores_files, arity: '3')

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    plot_factor_distributions.R F1_CCD3F.txt F2_CCD3F.txt F3_CCD3F.txt
    """

    stub:
    """
    touch distribution_diagnostics.log
    touch factor_distribution_F1.png
    touch factor_distribution_F2.png
    touch factor_distribution_F3.png
    """

    output:
    path("distribution_diagnostics.log"),   emit: log
    path("factor_distribution_F1.png"),     emit: f1_plot
    path("factor_distribution_F2.png"),     emit: f2_plot
    path("factor_distribution_F3.png"),     emit: f3_plot
}
