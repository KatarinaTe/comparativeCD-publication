// ============================================================================
// main.nf — 4. Plotting workflow
//
// Terminal, diagnostic-only plots over 3_Merging_filtering's outputs: PCA colored by breed,
// factor-score distributions, and per-item/per-factor response histograms. Nothing here feeds
// any further stage.
//
// Source: 02_data_processing/4_Plotting/{plot_PCA.R, plot_factor_distr.R,
//          plot_item&factor_hist.R}
// ============================================================================

include { validateParameters } from 'plugin/nf-schema'

include { PLOT_PCA              } from './process_definitions.nf'
include { PLOT_FACTOR_DISTR     } from './process_definitions.nf'
include { PLOT_ITEM_FACTOR_HIST } from './process_definitions.nf'

// The 14 dogCD survey items, in the fixed order PLOT_ITEM_FACTOR_HIST's script expects —
// matches 3_Merging_filtering's own getItemIds().
def getItemIds() {
    ['7', '93', '95', '145', '146', '147', '148', '149', '150', '151', '152', '153', '154', '155']
}

workflow {

    main:
    validateParameters()

    ch_eigenvec       = channel.value(file(params.eigenvec,       checkIfExists: true))
    ch_eigenval       = channel.value(file(params.eigenval,       checkIfExists: true))
    ch_data6          = channel.value(file(params.data6,          checkIfExists: true))
    ch_dog_breeds_csv = channel.value(file(params.dog_breeds_csv, checkIfExists: true))

    pca_out = PLOT_PCA(ch_eigenvec, ch_eigenval, ch_data6, ch_dog_breeds_csv)

    // 3_Merging_filtering's published QC6 fam files — same naming convention as that pipeline's
    // own FILTER_FACTOR_QC6/FILTER_ITEM_QC6 outputs, referenced by explicit filename here rather
    // than a glob to keep item order guaranteed and validated up front (checkIfExists: true).
    ch_f1_fam = channel.value(file("${params.qc6_fam_dir}/CCDF1_DA_MERGED_GENCOVE_AXIOM_QC6.fam", checkIfExists: true))
    ch_f2_fam = channel.value(file("${params.qc6_fam_dir}/CCDF2_DA_MERGED_GENCOVE_AXIOM_QC6.fam", checkIfExists: true))
    ch_f3_fam = channel.value(file("${params.qc6_fam_dir}/CCDF3_DA_MERGED_GENCOVE_AXIOM_QC6.fam", checkIfExists: true))

    distr_out = PLOT_FACTOR_DISTR(ch_f1_fam, ch_f2_fam, ch_f3_fam)

    ch_item_fams = channel.value(
        getItemIds().collect { id -> file("${params.qc6_fam_dir}/CCDitem${id}_DA_MERGED_GENCOVE_AXIOM_QC6.fam", checkIfExists: true) }
    )

    hist_out = PLOT_ITEM_FACTOR_HIST(ch_f1_fam, ch_f2_fam, ch_f3_fam, ch_item_fams)

    publish:
    pca_log               = pca_out.log
    pca_by_breed_plot     = pca_out.pca_by_breed
    factor_distr_log      = distr_out.log
    factor_distr_plot     = distr_out.plot
    item_factor_hist_log  = hist_out.log
    factor_hist_row_plot  = hist_out.factor_hist_row
    item_hist_row_f1_plot = hist_out.item_hist_row_f1
    item_hist_row_f2_plot = hist_out.item_hist_row_f2
    item_hist_row_f3_plot = hist_out.item_hist_row_f3

}

// Output publishing mode is symlink by default. This saves space over copying.
output {
    pca_log {
        path 'diagnostics'
    }
    pca_by_breed_plot {
        path 'plots'
    }
    factor_distr_log {
        path 'diagnostics'
    }
    factor_distr_plot {
        path 'plots'
    }
    item_factor_hist_log {
        path 'diagnostics'
    }
    factor_hist_row_plot {
        path 'plots'
    }
    item_hist_row_f1_plot {
        path 'plots'
    }
    item_hist_row_f2_plot {
        path 'plots'
    }
    item_hist_row_f3_plot {
        path 'plots'
    }
}
