// ============================================================================
// process_definitions.nf — 06. cCRE — cre_overlap (Fig. 4e; EpicDog tissue specificity for the Results text)
//
// Source: 06_cCRE/testing_CRE_overlap_260617.R -- a 1179-line, largely exploratory script with
// many superseded duplicate sections. Only the two analyses the manuscript reports are reproduced
// here: the UU per-region odds ratios (Fig. 4e; was Fig. 4f until Figure4_261006 dropped the
// tissue-specificity panel) and the EpicDog tissue-specificity odds ratios (now Results text and
// Supplementary Table 13 only). The source script's own multi-panel rates/size figures aren't
// part of the manuscript. See bin/*.R headers for the full provenance/conversion story of each.
//
// BUILD_EPICDOG_TISSUE_OVERLAPS and BUILD_UU_REGION_BP_OVERLAPS both reuse the bedtools+base-R
// container already frozen for BUILD_CLUMP_REGIONS_100KB (env/06_cCRE/clump_regions.yml) --
// nothing here needs anything beyond bedtools + awk. PLOT_TISSUE_SPECIFICITY and
// PLOT_REGION_FOREST reuse EMISSION_STATE_PLOT's container (ggplot2/dplyr already pinned there;
// no tidyr/patchwork/scales needed for these single-panel forest plots).
// ============================================================================

process BUILD_EPICDOG_TISSUE_OVERLAPS {
    label 'process_single'
    container 'community.wave.seqera.io/library/clump_regions:09876be1daa6bcb4'

    input:
    path(epicdog_states_dir, arity: '1')
    path(dogcd_gwas_regions, arity: '1')
    path(size_gwas_regions, arity: '1')

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    EPICDOG_DIR=${epicdog_states_dir} DOGCD_GWAS=${dogcd_gwas_regions} SIZE_GWAS=${size_gwas_regions} \\
        OUTPUT=epicdog_tissue_bp_overlaps.txt OUTDIR=bedtools_tmp \\
        check_epicdog_bp_overlap.sh
    """

    stub:
    """
    touch epicdog_tissue_bp_overlaps.txt
    """

    output:
    path("epicdog_tissue_bp_overlaps.txt"), emit: overlaps
}


process PLOT_TISSUE_SPECIFICITY {
    // Odds ratio (Brain vs. other tissues), dogCD and dogSize, EpicDog resource. Reported in the
    // Results text and Supplementary Table 13; no longer a figure panel (dropped from Fig. 4, 2026-10-06).
    label 'process_single'
    container 'community.wave.seqera.io/library/emission_state:0737d6956e5a4a7d'

    input:
    path(epicdog_tissue_overlaps, arity: '1')

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    plot_tissue_specificity.R ${epicdog_tissue_overlaps} tissue_specificity.pdf
    """

    stub:
    """
    touch tissue_specificity.pdf
    """

    output:
    path("tissue_specificity.pdf"), emit: plot
}


process BUILD_UU_REGION_BP_OVERLAPS {
    label 'process_single'
    container 'community.wave.seqera.io/library/clump_regions:09876be1daa6bcb4'

    input:
    path(uu_states_dir, arity: '1')
    path(dogcd_gwas_regions, arity: '1')
    path(size_gwas_regions, arity: '1')

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    STATES_DIR=${uu_states_dir} DOGCD_GWAS=${dogcd_gwas_regions} SIZE_GWAS=${size_gwas_regions} \\
        DOGCD_OUTPUT=brain_region_bp_overlaps.txt SIZE_OUTPUT=SIZE_brain_region_bp_overlaps.txt \\
        OUTDIR=bedtools_tmp \\
        check_uu_region_bp_overlap.sh
    """

    stub:
    """
    touch brain_region_bp_overlaps.txt SIZE_brain_region_bp_overlaps.txt
    """

    output:
    path("brain_region_bp_overlaps.txt"), emit: dogcd_overlaps
    path("SIZE_brain_region_bp_overlaps.txt"), emit: size_overlaps
}


process PLOT_REGION_FOREST {
    // Fig. 4e: odds ratio (dogCD vs. dogSize) per UU brain region (x8), ACG highlighted.
    label 'process_single'
    container 'community.wave.seqera.io/library/emission_state:0737d6956e5a4a7d'

    input:
    path(dogcd_region_overlaps, arity: '1')
    path(size_region_overlaps, arity: '1')
    path(dogcd_gwas_regions, arity: '1')
    path(size_gwas_regions, arity: '1')

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    dogcd_bp=\$(grep -v '^#' ${dogcd_gwas_regions} | awk '{sum+=\$3-\$2} END{print sum}')
    size_bp=\$(awk '{sum+=\$3-\$2} END{print sum}' ${size_gwas_regions})
    plot_region_forest.R ${dogcd_region_overlaps} ${size_region_overlaps} \${dogcd_bp} \${size_bp} \\
        region_forest.pdf region_OR_summary.txt
    """

    stub:
    """
    touch region_forest.pdf region_OR_summary.txt
    """

    output:
    path("region_forest.pdf"), emit: plot
    path("region_OR_summary.txt"), emit: summary
}
