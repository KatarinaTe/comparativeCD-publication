// ============================================================================
// main.nf — 04. GWAS — 5. Plotting
//
// Manhattan + QQ plots for Fig 1c: all 4 mlma-loco phenotypes (CCDF1-3, SIZE) and all 14
// polmm survey items.
//
// Source: 07_additional_plots_and_analyses/plot_gwas_results.R (moved here, see
// process_definitions.nf header)
// ============================================================================

include { validateParameters } from 'plugin/nf-schema'

include { PLOT_GWAS_RESULTS } from './process_definitions.nf'

def getMlmaPhenoIds()  { ['CCDF1', 'CCDF2', 'CCDF3', 'SIZE'] }
def getPolmmItemIds()  { ['7', '93', '95', '145', '146', '147', '148', '149', '150', '151', '152', '153', '154', '155'] }

workflow {

    main:
    validateParameters()

    ch_mlma_files = channel.fromList(getMlmaPhenoIds())
        .map { pheno -> file("${params.mlma_dir}/${pheno}_*_QC6_LOCO.loco.mlma", checkIfExists: true) }
        .collect()

    ch_polmm_files = channel.fromList(getPolmmItemIds())
        .map { id -> file("${params.polmm_dir}/modi_simuMarkerOutput_POLMM_item${id}_FULLGRM_fullGeno.txt", checkIfExists: true) }
        .collect()

    plot_out = PLOT_GWAS_RESULTS(ch_mlma_files, ch_polmm_files)

    publish:
    manhattan = plot_out.manhattan
    qq        = plot_out.qq

}

// Output publishing mode is symlink by default. This saves space over copying.
output {
    manhattan {
        path 'plotting'
    }
    qq {
        path 'plotting'
    }
}
