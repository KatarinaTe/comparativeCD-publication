// ============================================================================
// process_definitions.nf — 04. GWAS — 1. MLMA-LOCO (factors + SIZE)
//
// Source: 04_gwas/1_mlma-loco/{MLMA_PER_FACTOR.sh, clump_plink_factors.sh,
//          export_gwas_sumstats_CCD_MLMA.R}
//
// Genericized across F1/F2/F3 via a config channel, plus SIZE run through the identical
// mlma-loco step, clumping, and sumstats export. CLUMP_FACTOR/EXPORT_SUMSTATS branch their output
// prefix on pheno_id ('SIZE' vs. "CCD${pheno_id}") rather than hardcoding "CCD", matching
// 1_mlma-loco.nf's own bfilePrefix() convention for the same SIZE-has-no-"CCD"-prefix quirk.
//
// RUN_MLMA_LOCO's --out uses the full QC6 bfile prefix (e.g. CCDF1_DA_MERGED_GENCOVE_AXIOM_QC6_LOCO)
// — unambiguous across all 4 phenotypes in a flat published output directory; the filename itself
// is not parsed downstream.
//
// `--thread-num` uses `task.cpus` rather than a hardcoded value — thread count doesn't affect
// GWAS results, only wall-clock time.
//
// Reuses 3_Merging_filtering's container (gcta64 + R/dplyr) for mlma-loco/sumstats export, and the
// existing plink container for clumping — no new build needed.
// ============================================================================

process RUN_MLMA_LOCO {
    tag   "${factor_id}"
    label 'process_high'
    container 'community.wave.seqera.io/library/merging_filtering_tools:0e56b6d8b7fd7863'

    input:
    tuple val(factor_id), path(qc6_bed), path(qc6_bim), path(qc6_fam), path(phen), path(qcovar), path(covar)

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = qc6_bed.baseName
    """
    gcta64 --mlma-loco --bfile ${prefix} --pheno ${phen} --out ${prefix}_LOCO --thread-num ${task.cpus} --qcovar ${qcovar} --covar ${covar} --autosome-num 38
    """

    stub:
    """
    touch ${qc6_bed.baseName}_LOCO.loco.mlma
    touch ${qc6_bed.baseName}_LOCO.log
    """

    output:
    tuple val(factor_id), path("${qc6_bed.baseName}_LOCO.loco.mlma"), emit: mlma
    path("${qc6_bed.baseName}_LOCO.log"), emit: log
}


process CLUMP_FACTOR {
    // Runs for all 4 phenotypes (F1/F2/F3 + SIZE) — see this file's header.
    tag   "${pheno_id}"
    label 'process_medium'
    container 'community.wave.seqera.io/library/plink:1.90b7.7--feb905bcf3811778'

    input:
    tuple val(pheno_id), path(qc6_bed), path(qc6_bim), path(qc6_fam), path(mlma)
    path(gene_range_file, arity: '1')

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = qc6_bed.baseName
    def label  = pheno_id == 'SIZE' ? 'SIZE' : "CCD${pheno_id}"
    """
    plink --dog --bfile ${prefix} --out ${label}_clump250 --clump ${mlma} --clump-field p --clump-p1 0.000001 --clump-p2 0.00001 --clump-kb 250 --clump-r2 0.50 --clump-range ${gene_range_file}
    """

    stub:
    def label = pheno_id == 'SIZE' ? 'SIZE' : "CCD${pheno_id}"
    """
    touch ${label}_clump250.clumped
    touch ${label}_clump250.clumped.ranges
    touch ${label}_clump250.log
    """

    output:
    path("${pheno_id == 'SIZE' ? 'SIZE' : "CCD${pheno_id}"}_clump250.clumped"), emit: clumped
    path("${pheno_id == 'SIZE' ? 'SIZE' : "CCD${pheno_id}"}_clump250.clumped.ranges"), emit: clumped_ranges
    path("${pheno_id == 'SIZE' ? 'SIZE' : "CCD${pheno_id}"}_clump250.log"), emit: log
}


process EXPORT_SUMSTATS {
    // Runs for all 4 phenotypes (F1/F2/F3 + SIZE) — see this file's header.
    tag   "${pheno_id}"
    label 'process_single'
    container 'community.wave.seqera.io/library/merging_filtering_tools:0e56b6d8b7fd7863'

    input:
    tuple val(pheno_id), path(mlma), path(phen)

    when:
    task.ext.when == null || task.ext.when

    script:
    def label = pheno_id == 'SIZE' ? 'SIZE' : "CCD${pheno_id}"
    """
    export_gwas_sumstats_mlma.R ${mlma} ${phen} ${label}_sumstats.txt
    """

    stub:
    def label = pheno_id == 'SIZE' ? 'SIZE' : "CCD${pheno_id}"
    """
    touch ${label}_sumstats.txt
    """

    output:
    path("${pheno_id == 'SIZE' ? 'SIZE' : "CCD${pheno_id}"}_sumstats.txt"), emit: sumstats
}
