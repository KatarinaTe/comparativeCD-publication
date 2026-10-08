// ============================================================================
// process_definitions.nf — 04. GWAS — 2. POLMM (items)
//
// Source: 04_gwas/2_polmm/{POLMMgrab_PER_ITEM.sh, clump_plink_items.sh,
//          export_gwas_sumstats_CCD_POLMM.R}
//
// Genericized across all 14 items via a config channel.
//
// New container (community.wave.seqera.io/library/polmm_grab) — the GRAB R package (POLMM
// method) isn't part of 3_Merging_filtering's existing R stack. Clumping reuses the existing
// plink container; sumstats export's R script needs no packages beyond base R.
//
// GRAB is pinned to 0.2.5, the newest version available on conda-forge/CRAN; the exact version
// used in the original analysis is unconfirmed.
// ============================================================================

process RUN_POLMM_ITEM {
    // Null model fit (on the pruned EigenInput genotype) + full-genome marker test + the
    // original's own awk reformat of GRAB.Marker's colon-packed ID columns. See
    // bin/run_polmm_item.R and bin/reformat_polmm_output.sh.
    tag   "item${item_id}"
    label 'process_high'
    container 'community.wave.seqera.io/library/polmm_grab:4c906c4feeeb8621'

    input:
    tuple val(item_id), path(grab_file), path(qc6_bed), path(qc6_bim), path(qc6_fam), path(eigen_bed), path(eigen_bim), path(eigen_fam)

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    run_polmm_item.R ${grab_file} ${qc6_bed} ${eigen_bed} ${item_id} simuMarkerOutput_POLMM_item${item_id}_FULLGRM_fullGeno.txt
    reformat_polmm_output.sh simuMarkerOutput_POLMM_item${item_id}_FULLGRM_fullGeno.txt modi_simuMarkerOutput_POLMM_item${item_id}_FULLGRM_fullGeno.txt
    """

    stub:
    """
    touch simuMarkerOutput_POLMM_item${item_id}_FULLGRM_fullGeno.txt
    touch modi_simuMarkerOutput_POLMM_item${item_id}_FULLGRM_fullGeno.txt
    """

    output:
    tuple val(item_id), path("modi_simuMarkerOutput_POLMM_item${item_id}_FULLGRM_fullGeno.txt"), emit: modi_sumstats
    path("simuMarkerOutput_POLMM_item${item_id}_FULLGRM_fullGeno.txt"), emit: raw_marker
}


process CLUMP_ITEM {
    // All 14 items — unlike 1_mlma-loco's clumping, which only covered F1/F2/F3.
    tag   "item${item_id}"
    label 'process_medium'
    container 'community.wave.seqera.io/library/plink:1.90b7.7--feb905bcf3811778'

    input:
    tuple val(item_id), path(qc6_bed), path(qc6_bim), path(qc6_fam), path(modi_sumstats)
    path(gene_range_file, arity: '1')

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = qc6_bed.baseName
    """
    plink --dog --bfile ${prefix} --out polmmCCDitem${item_id}_clump250 --clump ${modi_sumstats} --clump-field Pvalue --clump-p1 0.000001 --clump-p2 0.00001 --clump-kb 250 --clump-r2 0.50 --clump-range ${gene_range_file}
    """

    stub:
    """
    touch polmmCCDitem${item_id}_clump250.clumped
    touch polmmCCDitem${item_id}_clump250.clumped.ranges
    touch polmmCCDitem${item_id}_clump250.log
    """

    output:
    path("polmmCCDitem${item_id}_clump250.clumped"), emit: clumped
    path("polmmCCDitem${item_id}_clump250.clumped.ranges"), emit: clumped_ranges
    path("polmmCCDitem${item_id}_clump250.log"), emit: log
}


process EXPORT_SUMSTATS_POLMM {
    tag   "item${item_id}"
    label 'process_single'
    container 'community.wave.seqera.io/library/merging_filtering_tools:0e56b6d8b7fd7863'

    input:
    tuple val(item_id), path(modi_sumstats)

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    export_gwas_sumstats_polmm.R ${modi_sumstats} ITEM${item_id}_sumstats.txt
    """

    stub:
    """
    touch ITEM${item_id}_sumstats.txt
    """

    output:
    path("ITEM${item_id}_sumstats.txt"), emit: sumstats
}
