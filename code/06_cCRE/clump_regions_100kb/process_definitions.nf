// ============================================================================
// process_definitions.nf — 06. cCRE — GWAS clump-region derivation
//
// Split out of the main cCRE.nf into its own substage (2026-09-15) so it can be run and
// stub-tested independently of EMISSION_STATE_PLOT — the two are unrelated analyses that only
// briefly shared one workflow file, and bundling them meant a missing 04_gwas/1_mlma-loco output
// broke EMISSION_STATE_PLOT's own, unrelated smoke test too.
//
// build_clump_regions_100kb.R (see bin/ header for the full method) — one BED region per
// genome-wide significant GWAS clump (04_gwas/1_mlma-loco's CLUMP_FACTOR output, plus
// 04_gwas/2_polmm's CLUMP_ITEM output for dogCD), spanning that clump's full extent (lead +
// every secondary SNP) padded 100kb each side, pooled across every phenotype passed in and
// merged where overlapping. Called once for "dogCD" (CCDF1+CCDF2+CCDF3 + all 14 POLMM items
// pooled) and once for "SIZE" alone — see clump_regions_100kb.nf.
//
// Provides the GWAS-region input ../cre_overlap/'s BUILD_EPICDOG_TISSUE_OVERLAPS and
// BUILD_UU_REGION_BP_OVERLAPS each take a copy of (params.dogcd_gwas_regions/size_gwas_regions).
// Empirically validated 2026-09-23 against real historical GWAS output — see this substage's
// README Notes and REPRODUCIBILITY_AUDIT.md for the full result (dogCD: exact match to the
// deposited reference; SIZE: matches to within a small, independently-explained residual).
//
// label bumped process_single -> process_medium (2026-09-23): the "dogCD" pool now reads up to
// 14 POLMM modi_simuMarkerOutput_POLMM_item*_FULLGRM_fullGeno.txt files, each ~1GB on disk with
// ~10M rows -- R's read.table() overhead on tables that size repeatedly OOM-killed this process
// (exit 137, even after 3 retries) at process_single's 2-6GB range. Found running this for real
// against historical POLMM data on the cluster.
// ============================================================================

process BUILD_CLUMP_REGIONS_100KB {
    tag   "${label}"
    label 'process_medium'
    container 'community.wave.seqera.io/library/clump_regions:09876be1daa6bcb4'

    input:
    tuple val(label), path(clumped_and_mlma_files)

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    build_clump_regions_100kb.R unsorted.bed ${clumped_and_mlma_files}
    bedtools sort -i unsorted.bed | bedtools merge -i - > ${label}_gws100kb_regions.bed
    """

    stub:
    """
    touch ${label}_gws100kb_regions.bed
    """

    output:
    tuple val(label), path("${label}_gws100kb_regions.bed"), emit: regions
}
