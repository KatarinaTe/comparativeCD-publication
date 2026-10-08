// ============================================================================
// process_definitions.nf — 04. GWAS — 5. Plotting
//
// Source: 07_additional_plots_and_analyses/plot_gwas_results.R, moved here — it plots
// 1_mlma-loco's and 2_polmm's own outputs, so it belongs alongside its real parent stage
// (same reasoning as the 2026-09-15 GO-semantic-clustering move into 05_network).
//
// The script already reads/writes by relative filename (Sys.glob by phenotype prefix for
// mlma, exact filename by item id for polmm) — Nextflow stages the input paths into the task
// work dir under their original names, so no CLI-arg parameterization was needed beyond the
// shebang/bin/ move itself.
// ============================================================================

process PLOT_GWAS_RESULTS {
    label 'process_single'
    container 'community.wave.seqera.io/library/gwas_supplement_plots:e626727dcb23673b'

    input:
    path(mlma_files)
    path(polmm_files)

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    plot_gwas_results.R
    """

    stub:
    """
    for pheno in CCDF1 CCDF2 CCDF3 SIZE; do
        touch "\${pheno}_loco_manhattan.jpg" "\${pheno}_loco_qq.jpg"
    done
    for item in 7 93 95 145 146 147 148 149 150 151 152 153 154 155; do
        touch "CCDitem\${item}_polmm_manhattan.jpg" "CCDitem\${item}_polmm_qq.jpg"
    done
    """

    output:
    path("*_manhattan.jpg"), emit: manhattan
    path("*_qq.jpg"),        emit: qq
}
