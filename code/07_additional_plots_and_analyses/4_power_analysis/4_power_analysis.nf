// ============================================================================
// main.nf — 07. Additional plots and analyses — 4. Power analysis
//
// Statistical power comparison, dogCD GWAS vs. human OCD GWAS (Strom et al. 2025 meta-analysis).
// DogCD achieves comparable power despite a far smaller sample size.
//
// Source: power_calcs.dog_v_human_GWAS.R (moved here, see process_definitions.nf header)
// ============================================================================

include { validateParameters } from 'plugin/nf-schema'

include { POWER_ANALYSIS } from './process_definitions.nf'

workflow {

    main:
    validateParameters()

    ch_dog_gwas_file   = channel.fromPath(params.dog_gwas_file, checkIfExists: true)
    ch_human_gwas_file = channel.fromPath(params.human_gwas_file, checkIfExists: true)

    power_out = POWER_ANALYSIS(ch_dog_gwas_file, ch_human_gwas_file)

    publish:
    dog_summary   = power_out.dog_summary
    human_summary = power_out.human_summary
    plot_png      = power_out.plot_png
    plot_pdf      = power_out.plot_pdf

}

// Output publishing mode is symlink by default. This saves space over copying.
output {
    dog_summary {
        path 'power_analysis'
    }
    human_summary {
        path 'power_analysis'
    }
    plot_png {
        path 'power_analysis'
    }
    plot_pdf {
        path 'power_analysis'
    }
}
