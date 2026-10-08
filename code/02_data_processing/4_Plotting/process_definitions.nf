// ============================================================================
// process_definitions.nf — 4. Plotting
//
// Source: 02_data_processing/4_Plotting/{plot_PCA.R, plot_factor_distr.R,
//          plot_item&factor_hist.R}
//
// Reuses 3_Merging_filtering's container (community.wave.seqera.io/library/
// merging_filtering_tools) — same R/dplyr/ggplot2/patchwork/ggtext stack, no new build needed.
//
// Every diagnostic print in the original scripts is reproduced as a real sink()-captured log.
// Several values baked into the *visible* plots in the original (PCA variance-% axis labels, N=
// annotations, axis ranges for continuous scores) are computed from the real data here instead.
// See each bin/*.R header for exactly what's computed vs. kept literal.
// ============================================================================

process PLOT_PCA {
    // PCA scatter plot colored/shaped by breed. Breed metadata comes from dog_breeds_csv
    // (data/02_data_processing/DarwinsArk_20220715_dogs_genotyped_breed_sex.csv) — see bin/plot_pca.R's header for
    // the variance-% axis-label computation and why the original's uncolored first plot (whose
    // legend didn't match its own data) is dropped rather than reproduced.
    tag   "cohort"
    label 'process_medium'
    container 'community.wave.seqera.io/library/merging_filtering_tools:0e56b6d8b7fd7863'

    input:
    path(eigenvec, arity: '1')
    path(eigenval, arity: '1')
    path(data6, arity: '1')
    path(dog_breeds_csv, arity: '1')

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    plot_pca.R ${eigenvec} ${eigenval} ${data6} ${dog_breeds_csv}
    """

    stub:
    """
    touch pca_diagnostics.log
    touch pca_by_breed.png
    """

    output:
    path("pca_diagnostics.log"), emit: log
    path("pca_by_breed.png"),    emit: pca_by_breed
}


process PLOT_FACTOR_DISTR {
    // F1/F2/F3 factor-score distribution histograms, stacked vertically. See
    // bin/plot_factor_distr.R's header for the per-factor config (only F3's block styles
    // axis.title.x) and the dynamically-computed N=/shared-axis-range values.
    tag   "cohort"
    label 'process_single'
    container 'community.wave.seqera.io/library/merging_filtering_tools:0e56b6d8b7fd7863'

    input:
    path(f1_fam, arity: '1')
    path(f2_fam, arity: '1')
    path(f3_fam, arity: '1')

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    plot_factor_distr.R ${f1_fam} ${f2_fam} ${f3_fam}
    """

    stub:
    """
    touch factor_distr_diagnostics.log
    touch factor_distr.png
    """

    output:
    path("factor_distr_diagnostics.log"), emit: log
    path("factor_distr.png"),             emit: plot
}


process PLOT_ITEM_FACTOR_HIST {
    // Factor-score histogram row (F1+F2+F3) and three per-item response histogram rows (CCDF1's
    // 4 items, CCDF2's 3 items, CCDF3's 7 items). See bin/plot_item_factor_hist.R's header for
    // the per-plot config tables (title, color, x-axis breaks/margins — reproduced exactly per
    // plot, including item152's missing scale_x_continuous() call, a genuine source oddity kept
    // as-is) and item147's `c(0.5,6,5)` -> `c(0.5,6.5)` fix.
    //
    // item_fams must be exactly the 14 items in this fixed order — matches the order
    // plot_item_factor_hist.R's arg-parsing expects, not alphabetical/numeric:
    // 7, 93, 95, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155.
    tag   "cohort"
    label 'process_single'
    container 'community.wave.seqera.io/library/merging_filtering_tools:0e56b6d8b7fd7863'

    input:
    path(f1_fam, arity: '1')
    path(f2_fam, arity: '1')
    path(f3_fam, arity: '1')
    path(item_fams, arity: '14')

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    plot_item_factor_hist.R ${f1_fam} ${f2_fam} ${f3_fam} ${item_fams.join(' ')}
    """

    stub:
    """
    touch item_factor_hist_diagnostics.log
    touch factor_hist_row.png
    touch item_hist_row_f1.png
    touch item_hist_row_f2.png
    touch item_hist_row_f3.png
    """

    output:
    path("item_factor_hist_diagnostics.log"), emit: log
    path("factor_hist_row.png"),              emit: factor_hist_row
    path("item_hist_row_f1.png"),             emit: item_hist_row_f1
    path("item_hist_row_f2.png"),             emit: item_hist_row_f2
    path("item_hist_row_f3.png"),             emit: item_hist_row_f3
}
