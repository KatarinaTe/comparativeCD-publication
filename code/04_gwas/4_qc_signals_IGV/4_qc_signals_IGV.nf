// ============================================================================
// main.nf — 04. GWAS — 4. Manual IGV signal QC — automated prep tooling
//
// Automates everything leading up to a human's IGV read-pileup check: extracts one lead SNP
// (lowest P) per finemapping credible set, maps BAM-having individuals to their genotype-ID
// (dogID), extracts their genotypes at those lead SNPs, and picks up to N candidate individuals
// per genotype class to actually load into IGV — re-runnable per --round for a fresh batch.
//
// Source: 04_gwas/4_qc_signals_IGV/check_signals_igv_manually.sh (informal scratch notes, not a
// finished script).
// ============================================================================

include { validateParameters } from 'plugin/nf-schema'

include { EXTRACT_CREDIBLE_SET_SNPS  } from './process_definitions.nf'
include { BUILD_BAM_DOGID_KEEP_LIST  } from './process_definitions.nf'
include { EXTRACT_LEAD_SNP_GENOTYPES } from './process_definitions.nf'
include { SELECT_IGV_INDIVIDUALS     } from './process_definitions.nf'

workflow {

    main:
    validateParameters()

    // ── Lead-SNP list: top-P SNP per credible set, across every 3_finemap_susie result file ──

    ch_finemap_results = channel.fromPath("${params.finemap_results_dir}/*.gz", checkIfExists: true)
        .collect()

    snps_out = EXTRACT_CREDIBLE_SET_SNPS(ch_finemap_results)

    // ── BAM-having sampleIDs -> dogID, via the already-deposited QC3modi crosswalk fam ──

    ch_bam_listing = channel.fromPath("${params.bam_dir}/*.bam", checkIfExists: true)
        .map { it.name }
        .collectFile(name: 'bam_listing.txt', newLine: true)

    ch_crosswalk_fam = channel.value(file(params.crosswalk_fam, checkIfExists: true))

    keep_out = BUILD_BAM_DOGID_KEEP_LIST(ch_bam_listing, ch_crosswalk_fam)

    // ── Genotype extraction at the lead SNPs, restricted to BAM-having dogs ──

    ch_qc5_bed = channel.value(file("${params.qc6_plink_dir}/DA_MERGED_GENCOVE_AXIOM_QC5.bed", checkIfExists: true))
    ch_qc5_bim = channel.value(file("${params.qc6_plink_dir}/DA_MERGED_GENCOVE_AXIOM_QC5.bim", checkIfExists: true))
    ch_qc5_fam = channel.value(file("${params.qc6_plink_dir}/DA_MERGED_GENCOVE_AXIOM_QC5.fam", checkIfExists: true))
    ch_qc5_plink = ch_qc5_bed.combine(ch_qc5_bim).combine(ch_qc5_fam)

    geno_out = EXTRACT_LEAD_SNP_GENOTYPES(ch_qc5_plink, snps_out.snp_list, keep_out.keep_list)

    // ── Per-genotype-class individual selection — re-run with --round N for a fresh batch ──

    selection_out = SELECT_IGV_INDIVIDUALS(geno_out.raw, params.round, params.n_per_class, params.seed)

    publish:
    lead_snps           = snps_out.snp_list
    bam_dogid_keep_list = keep_out.keep_list
    unmatched_bam_ids   = keep_out.unmatched
    lead_snp_genotypes  = geno_out.raw
    lead_snp_geno_log   = geno_out.log
    igv_selection       = selection_out.selection

}

// Output publishing mode is symlink by default. This saves space over copying.
output {
    lead_snps {
        path 'lead_snps'
    }
    bam_dogid_keep_list {
        path 'bam_mapping'
    }
    unmatched_bam_ids {
        path 'bam_mapping'
    }
    lead_snp_genotypes {
        path 'genotypes'
    }
    lead_snp_geno_log {
        path 'genotypes'
    }
    igv_selection {
        path 'selection'
    }
}
