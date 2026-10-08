include { CONCAT_CHR_VCFS         } from '../process_definitions.nf'
include { GLIMPSE_LIKELIHOODS_DOWN   } from '../process_definitions.nf'
include { GLIMPSE_IMPUTE_CHUNKS_DOWN } from '../process_definitions.nf'
include { INDEX_BCFS_DOWN            } from '../process_definitions.nf'
include { GLIMPSE_LIGATE_DOWN        } from '../process_definitions.nf'
include { BCF_VCF_DOWN               } from '../process_definitions.nf'
include { GLIMPSE_QC_FILTER_DOWN     } from '../process_definitions.nf'

// Re-impute every downsampling fraction via GLIMPSE, in one call. The *_DOWN processes are
// forks of 01_mapping's GLIMPSE processes with a `meta` (fraction) value threaded through
// every tuple, keeping different fractions' same-chromosome data apart through the chr-keyed
// grouping steps below.
workflow IMPUTE_DOWNSAMPLED_FRACTIONS {
    take:
    ch_bam_bundles_by_fraction // channel: [fraction, [bam...], [bai...]]
    chunks                     // channel: [chr, chunks_file] — fraction-independent
    sites                      // channel: [chr, sites_vcf, sites_csi, sites_tsv, sites_tbi] — fraction-independent
    ch_assembly_ref
    ch_glimpse_ref

    main:
    ch_for_gl = sites.combine(ch_bam_bundles_by_fraction)
        .map { chr, sites_vcf, sites_csi, sites_tsv, sites_tbi, fraction, bams, bais ->
            [fraction, chr, sites_vcf, sites_csi, sites_tsv, sites_tbi, bams, bais]
        }
    gl_out = GLIMPSE_LIKELIHOODS_DOWN(ch_for_gl, ch_assembly_ref)

    ch_chunk_records = chunks.flatMap { _chr, chunks_file ->
        chunks_file.readLines().collect { line ->
            def fields = line.trim().split(/\s+/)
            // GLIMPSE chunk format: ID CHR INPUT_REGION OUTPUT_REGION
            [fields[1], fields[0], fields[2], fields[3]]
        }
    }
    // Every chunk is imputed once per fraction.
    ch_gl_by_chr = gl_out.gl_vcf.map { fraction, chr, gl_vcf, gl_vcf_csi -> [chr, fraction, gl_vcf, gl_vcf_csi] }
    ch_for_impute = ch_chunk_records
        .combine(ch_gl_by_chr, by: 0)
        .map { chr, chunk_id, input_region, output_region, fraction, gl_vcf, gl_vcf_csi ->
            [fraction, chr, chunk_id, input_region, output_region, gl_vcf, gl_vcf_csi]
        }

    impute_out  = GLIMPSE_IMPUTE_CHUNKS_DOWN(ch_for_impute, ch_glimpse_ref)
    indexed_out = INDEX_BCFS_DOWN(impute_out.imputed_bcf)

    ch_for_ligate = indexed_out.imputed
        .map { fraction, chr, _chunk_id, bcf, csi -> [fraction, chr, bcf, csi] }
        .groupTuple(by: [0, 1])
    ligate_out = GLIMPSE_LIGATE_DOWN(ch_for_ligate)

    bcf_vcf_out = BCF_VCF_DOWN(ligate_out.merged_bcf.groupTuple(by: [0, 1]))
    // GLIMPSE_QC_FILTER_DOWN still runs here, matching 2_mapping&imputing&depth_bams.sh's Step 5 —
    // but 3_genotype_concordance.sh's own concat command reads the pre-QC-filter merged VCF
    // (bcf_vcf_out, below), not this step's qc_modi output, despite the output filename implying
    // otherwise. Preserved as-is rather than "corrected" — see the stage README.
    GLIMPSE_QC_FILTER_DOWN(bcf_vcf_out.merged)

    ch_for_concat = bcf_vcf_out.merged
        .map { fraction, _chr, vcf, _csi -> ["downsampled_${fraction}", vcf] }
        .groupTuple()

    concat_out = CONCAT_CHR_VCFS(ch_for_concat)

    emit:
    vcf = concat_out.vcf // [label, vcf, csi]
}
