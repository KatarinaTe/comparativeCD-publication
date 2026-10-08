// ============================================================================
// process_definitions.nf — 07. Additional plots and analyses — 4. Power analysis
//
// Source: power_calcs.dog_v_human_GWAS.R, moved here; the two hardcoded relative filenames
// (Strom.Table1.tsv, ST5.dogCD.top_gwas.tsv) became CLI args, no other logic changed. See the
// script's own header for the full provenance story.
// ============================================================================

process POWER_ANALYSIS {
    label 'process_single'
    container 'community.wave.seqera.io/library/go_semantic_plots:69629d5c597a400b'

    input:
    path(dog_gwas_file, arity: '1')
    path(human_gwas_file, arity: '1')

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    power_calcs.dog_v_human_GWAS.R ${dog_gwas_file} ${human_gwas_file}
    """

    stub:
    """
    touch dogCD_power_summary.csv
    touch human_OCD_power_summary.csv
    touch dogCD_vs_humanOCD.GWAS_power_comparison.png
    touch dogCD_vs_humanOCD.GWAS_power_comparison.pdf
    """

    output:
    path("dogCD_power_summary.csv"),                         emit: dog_summary
    path("human_OCD_power_summary.csv"),                      emit: human_summary
    path("dogCD_vs_humanOCD.GWAS_power_comparison.png"), emit: plot_png
    path("dogCD_vs_humanOCD.GWAS_power_comparison.pdf"), emit: plot_pdf
}
