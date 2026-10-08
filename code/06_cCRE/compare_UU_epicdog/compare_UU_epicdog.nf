// ============================================================================
// compare_UU_epicdog.nf — 06. cCRE — UU vs. EpicDog element-count comparison (Fig. 4d)
//
// Source: 06_cCRE/compare_UU_epicdog_elements.R, 06_cCRE/build_all_tissues_filtered.sh
// ============================================================================

include { validateParameters } from 'plugin/nf-schema'

include { BUILD_ALL_TISSUES_FILTERED } from './process_definitions.nf'
include { COMPARE_UU_EPICDOG_ELEMENTS } from './process_definitions.nf'

workflow {

    main:
    validateParameters()

    ch_epic_cerebellum_bed = channel.fromPath(params.epic_cerebellum_bed, checkIfExists: true)
    ch_epic_cerebrum_bed   = channel.fromPath(params.epic_cerebrum_bed, checkIfExists: true)

    ch_uu_acg_bed              = channel.fromPath(params.uu_acg_bed, checkIfExists: true)
    ch_uu_cerebellum_bed       = channel.fromPath(params.uu_cerebellum_bed, checkIfExists: true)
    ch_uu_frontal_lobe_bed     = channel.fromPath(params.uu_frontal_lobe_bed, checkIfExists: true)
    ch_uu_hypothalamus_bed     = channel.fromPath(params.uu_hypothalamus_bed, checkIfExists: true)
    ch_uu_occipital_cortex_bed = channel.fromPath(params.uu_occipital_cortex_bed, checkIfExists: true)
    ch_uu_striatum_bed         = channel.fromPath(params.uu_striatum_bed, checkIfExists: true)
    ch_uu_temporal_cortex_bed  = channel.fromPath(params.uu_temporal_cortex_bed, checkIfExists: true)
    ch_uu_thalamus_bed         = channel.fromPath(params.uu_thalamus_bed, checkIfExists: true)

    build_out = BUILD_ALL_TISSUES_FILTERED(
        ch_uu_acg_bed,
        ch_uu_cerebellum_bed,
        ch_uu_frontal_lobe_bed,
        ch_uu_hypothalamus_bed,
        ch_uu_occipital_cortex_bed,
        ch_uu_striatum_bed,
        ch_uu_temporal_cortex_bed,
        ch_uu_thalamus_bed
    )

    compare_out = COMPARE_UU_EPICDOG_ELEMENTS(ch_epic_cerebellum_bed, ch_epic_cerebrum_bed, build_out.bed)

    publish:
    all_tissues_bed = build_out.bed
    table           = compare_out.table
    pdf             = compare_out.pdf

}

// Output publishing mode is symlink by default. This saves space over copying.
output {
    all_tissues_bed {
        path 'compare_UU_epicdog'
    }
    table {
        path 'compare_UU_epicdog'
    }
    pdf {
        path 'compare_UU_epicdog'
    }
}
