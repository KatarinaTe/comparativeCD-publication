// ============================================================================
// clump_regions_100kb.nf — 06. cCRE — GWAS clump-region derivation workflow
//
// Split out of cCRE.nf (2026-09-15) into its own independent substage — see
// process_definitions.nf's header for why.
// ============================================================================

include { validateParameters } from 'plugin/nf-schema'

include { BUILD_CLUMP_REGIONS_100KB } from './process_definitions.nf'

workflow {

    main:
    validateParameters()

    // dogCD = CCDF1+CCDF2+CCDF3 (continuous factors, mlma-loco) plus all 14 POLMM item traits
    // pooled together; SIZE = its own single phenotype. Each entry is
    // [clumped_F1, mlma_F1, clumped_F2, mlma_F2, ...] — order matters, matches how
    // build_clump_regions_100kb.R reads its args in (clumped, mlma) pairs.
    def dogcd_ids = ['F1', 'F2', 'F3']
    def dogcd_factor_files = dogcd_ids.collectMany { id ->
        [
            file("${params.clumping_dir}/CCD${id}_clump250.clumped", checkIfExists: true),
            file("${params.mlma_dir}/CCD${id}_DA_MERGED_GENCOVE_AXIOM_QC6_LOCO.loco.mlma", checkIfExists: true)
        ]
    }
    // The 14 dogCD survey items — matches 2_polmm's own getItemIds().
    def item_ids = ['7', '93', '95', '145', '146', '147', '148', '149', '150', '151', '152', '153', '154', '155']
    def dogcd_item_files = item_ids.collectMany { id ->
        [
            file("${params.polmm_clumping_dir}/polmmCCDitem${id}_clump250.clumped", checkIfExists: true),
            file("${params.polmm_gwas_dir}/modi_simuMarkerOutput_POLMM_item${id}_FULLGRM_fullGeno.txt", checkIfExists: true)
        ]
    }
    def dogcd_files = dogcd_factor_files + dogcd_item_files
    def size_files = [
        file("${params.clumping_size_dir}/SIZE_clump250.clumped", checkIfExists: true),
        file("${params.mlma_dir}/SIZE_DA_MERGED_GENCOVE_AXIOM_QC6_LOCO.loco.mlma", checkIfExists: true)
    ]

    ch_region_input = channel.of(['dogCD', dogcd_files], ['SIZE', size_files])
    region_out = BUILD_CLUMP_REGIONS_100KB(ch_region_input)

    publish:
    regions = region_out.regions.map { _label, bed -> bed }

}

// Output publishing mode is symlink by default. This saves space over copying.
output {
    regions {
        path 'gws100kb_regions'
    }
}
