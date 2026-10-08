// ============================================================================
// process_definitions.nf — 07. Additional plots and analyses — 2. GWAS-catalog overlap
//
// Source: reproduce_gwas_catalog_overlap.R, moved here; hardcoded relative filenames replaced
// with CLI args, no other logic changed. See the script's own header for the full provenance
// story (a frozen stand-in for a Google-Sheet-only "gene list" tab with no file equivalent).
// ============================================================================

process REPRODUCE_GWAS_CATALOG_OVERLAP {
    label 'process_single'
    container 'community.wave.seqera.io/library/gwas_supplement_plots:e626727dcb23673b'

    input:
    path(catalog_hits_file, arity: '1')
    path(name_conversion_file, arity: '1')
    path(frozen_overlap_file, arity: '1')
    path(frozen_gene_list_file, arity: '1')

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    reproduce_gwas_catalog_overlap.R ${catalog_hits_file} ${name_conversion_file} ${frozen_overlap_file} ${frozen_gene_list_file}
    """

    stub:
    """
    touch signif_hits_in_gwas_catalog.csv
    touch gene_catalog_hit_summary.csv
    touch gene_list_with_catalog_hit_counts.csv
    touch gene_list_overlap_crosscheck.csv
    """

    output:
    path("signif_hits_in_gwas_catalog.csv"),       emit: signif_hits
    path("gene_catalog_hit_summary.csv"),          emit: hit_summary
    path("gene_list_with_catalog_hit_counts.csv"), emit: gene_list_with_hits
    path("gene_list_overlap_crosscheck.csv"),      emit: crosscheck
}
