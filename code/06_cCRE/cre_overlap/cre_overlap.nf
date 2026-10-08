// ============================================================================
// main.nf — 06. cCRE — cre_overlap workflow (Fig. 4e)
//
// Two independent analyses sharing GWAS-region inputs: EpicDog tissue-specificity (Results text and
// Supplementary Table 13; no longer a figure panel) and UU per-brain-region overlap (Fig. 4e). See process_definitions.nf's header for detail.
//
// Source: 06_cCRE/testing_CRE_overlap_260617.R
// ============================================================================

include { validateParameters } from 'plugin/nf-schema'

include { BUILD_EPICDOG_TISSUE_OVERLAPS } from './process_definitions.nf'
include { PLOT_TISSUE_SPECIFICITY       } from './process_definitions.nf'
include { BUILD_UU_REGION_BP_OVERLAPS   } from './process_definitions.nf'
include { PLOT_REGION_FOREST            } from './process_definitions.nf'

workflow {

    main:
    validateParameters()

    ch_epicdog_states_dir = channel.fromPath(params.epicdog_states_dir, checkIfExists: true)
    ch_uu_states_dir      = channel.fromPath(params.uu_states_dir, checkIfExists: true)
    ch_dogcd_gwas_regions = channel.fromPath(params.dogcd_gwas_regions, checkIfExists: true)
    ch_size_gwas_regions  = channel.fromPath(params.size_gwas_regions, checkIfExists: true)

    epicdog_out = BUILD_EPICDOG_TISSUE_OVERLAPS(ch_epicdog_states_dir, ch_dogcd_gwas_regions, ch_size_gwas_regions)
    tissue_plot_out = PLOT_TISSUE_SPECIFICITY(epicdog_out.overlaps)

    uu_out = BUILD_UU_REGION_BP_OVERLAPS(ch_uu_states_dir, ch_dogcd_gwas_regions, ch_size_gwas_regions)
    forest_plot_out = PLOT_REGION_FOREST(uu_out.dogcd_overlaps, uu_out.size_overlaps, ch_dogcd_gwas_regions, ch_size_gwas_regions)

    publish:
    epicdog_tissue_overlaps  = epicdog_out.overlaps
    tissue_specificity_plot  = tissue_plot_out.plot
    uu_dogcd_region_overlaps = uu_out.dogcd_overlaps
    uu_size_region_overlaps  = uu_out.size_overlaps
    region_forest_plot       = forest_plot_out.plot
    region_forest_summary    = forest_plot_out.summary

}

// Output publishing mode is symlink by default. This saves space over copying.
output {
    epicdog_tissue_overlaps {
        path 'overlaps'
    }
    tissue_specificity_plot {
        path 'tissue_specificity'
    }
    uu_dogcd_region_overlaps {
        path 'overlaps'
    }
    uu_size_region_overlaps {
        path 'overlaps'
    }
    region_forest_plot {
        path 'region_forest'
    }
    region_forest_summary {
        path 'region_forest'
    }
}
