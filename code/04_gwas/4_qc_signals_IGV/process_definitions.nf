// ============================================================================
// process_definitions.nf — 04. GWAS — 4. Manual IGV signal QC — automated prep tooling
//
// Source: 04_gwas/4_qc_signals_IGV/check_signals_igv_manually.sh — informal scratch notes, not a
// finished script to reproduce literally. This automates everything leading up to the one
// genuinely manual step (a human visually judging IGV read pileups): which SNP to check per
// credible set, which BAM-having individuals to load, and generating a fresh, non-overlapping
// batch of individuals if an earlier batch had no usable reads at that exact position.
//
// The lead-SNP list doesn't need any per-SNP cCRE/chromatin-state analysis at all: PolyFun's
// finemapper.py (SuSiE) already tags every SNP with a CREDIBLE_SET column in 3_finemap_susie's own
// finemap_results output (0 = not in a set).
//
// 01_mapping's BAM filenames use the raw sequencing sampleID; 3_Merging_filtering's QC5/QC6
// genotypes use "dogID" (relabeled via a positional fam substitution with no explicit lookup
// table published anywhere). DA_MERGED_GENCOVE_AXIOM_QC3modi.fam — already-deposited raw data,
// predating the relabel — carries both the original sampleID and the final dogID (IID) columns in
// the same row, and is reused directly as the crosswalk (see bin/build_bam_dogid_keep_list.py's
// header for the full reasoning).
//
// Reuses 3_finemap_susie's own polyfun container for all the Python steps (pandas already needed
// for the credible-set/genotype-selection logic; build_bam_dogid_keep_list.py only needs the
// stdlib, but reusing the same already-built container avoids sourcing a second one for one small
// script). Reuses the existing plink container for genotype extraction.
// ============================================================================

process EXTRACT_CREDIBLE_SET_SNPS {
    // See bin/extract_credible_set_snps.py — one SNP (lowest P) per credible set, across every
    // finemap_results/*.gz file from 3_finemap_susie.
    label 'process_single'
    container 'community.wave.seqera.io/library/polyfun:cb5ec2572b0e5a72'

    input:
    path(finemap_results)

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    extract_credible_set_snps.py ${finemap_results} lead_snps.txt
    """

    stub:
    """
    touch lead_snps.txt
    """

    output:
    path("lead_snps.txt"), emit: snp_list
}


process BUILD_BAM_DOGID_KEEP_LIST {
    // See bin/build_bam_dogid_keep_list.py's header for the sampleID<->dogID crosswalk reasoning.
    label 'process_single'
    container 'community.wave.seqera.io/library/polyfun:cb5ec2572b0e5a72'

    input:
    path(bam_listing, arity: '1')
    path(crosswalk_fam, arity: '1')

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    build_bam_dogid_keep_list.py ${bam_listing} ${crosswalk_fam} dogid_bam_keep.txt unmatched_bam_samples.txt
    """

    stub:
    """
    touch dogid_bam_keep.txt
    touch unmatched_bam_samples.txt
    """

    output:
    path("dogid_bam_keep.txt"), emit: keep_list
    path("unmatched_bam_samples.txt"), emit: unmatched
}


process EXTRACT_LEAD_SNP_GENOTYPES {
    // Whole-genome QC5 (not a per-chromosome split like 3_finemap_susie's shared genofile) —
    // simpler here since the lead SNPs span potentially all 38 chromosomes and there's no SuSiE
    // single-chromosome-scoping requirement driving a split the way there was in 3_finemap_susie.
    label 'process_medium'
    container 'community.wave.seqera.io/library/plink:1.90b7.7--feb905bcf3811778'

    input:
    tuple path(qc5_bed), path(qc5_bim), path(qc5_fam)
    path(lead_snps, arity: '1')
    path(keep_list, arity: '1')

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = qc5_bed.baseName
    """
    plink --dog --bfile ${prefix} --extract ${lead_snps} --keep ${keep_list} --recode A --out lead_snp_genotypes
    """

    stub:
    """
    touch lead_snp_genotypes.raw
    touch lead_snp_genotypes.log
    """

    output:
    path("lead_snp_genotypes.raw"), emit: raw
    path("lead_snp_genotypes.log"), emit: log
}


process SELECT_IGV_INDIVIDUALS {
    // See bin/select_igv_individuals.py — up to n_per_class individuals per genotype class per
    // lead SNP; re-run with a higher `round` (e.g. via CLI_OPTS='--round 2') for a fresh,
    // non-overlapping batch if an earlier round had no usable reads.
    label 'process_single'
    container 'community.wave.seqera.io/library/polyfun:cb5ec2572b0e5a72'

    input:
    path(raw_file, arity: '1')
    val(round)
    val(n_per_class)
    val(seed)

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    select_igv_individuals.py ${raw_file} igv_individual_selection.tsv --n-per-class ${n_per_class} --seed ${seed} --round ${round}
    """

    stub:
    """
    touch igv_individual_selection.tsv
    """

    output:
    path("igv_individual_selection.tsv"), emit: selection
}
