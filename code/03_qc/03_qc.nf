// ============================================================================
// main.nf — 03. QC workflow
//
// Three independent, one-off validation studies (none feed 04_gwas/05_network):
//   1. HighPass vs LowPass concordance   (11 dogs sequenced both ways)
//   2. Axiom-array vs LowPass concordance (5 duplicate-genotyped dogs + 1 sanity check)
//   3. Downsampling study                (10 HighPass dogs re-imputed at 10 target coverages)
//
// All three share one comparison tail (CONCORDANCE_ANALYSIS), called once over a channel
// mixing all 12 comparisons (study 1, study 2, 10 downsampling fractions) together, tagged by
// label. The downsampling study's GLIMPSE re-imputation chain uses forked copies of the
// 01_mapping processes (MERGE_BQSR_GVCF_DOWN, GLIMPSE_*_DOWN) with a `meta` fraction value
// threaded through every tuple, so it likewise runs once for all 10 fractions.
// ============================================================================

include { validateParameters; samplesheetToList } from 'plugin/nf-schema'

include { FILTER_PASS_VCF      } from './process_definitions.nf'
include { CONCAT_CHR_VCFS      } from './process_definitions.nf'
include { MAF_FILTER_VCF       } from './process_definitions.nf'
include { EXTRACT_SAMPLES_VCF  } from './process_definitions.nf'
include { DOWNSAMPLE_BAM       } from './process_definitions.nf'
include { MERGE_BQSR_GVCF_DOWN } from './process_definitions.nf'

include { GLIMPSE_DEFINE_CHUNKS  } from '../01_mapping/2_cohort_genotyping/process_definitions.nf'
include { GLIMPSE_VARIABLE_SITES } from '../01_mapping/2_cohort_genotyping/process_definitions.nf'

include { CONCORDANCE_ANALYSIS } from './subworkflows/concordance_analysis.nf'
include { IMPUTE_DOWNSAMPLED_FRACTIONS } from './subworkflows/impute_downsampled_fraction.nf'

workflow {

    main:
    validateParameters()

    def ref      = file(params.assembly_ref,                        checkIfExists: true)
    def ref_fai  = file(ref.resolveSibling("${ref.name}.fai"),      checkIfExists: true)
    def ref_dict = file(ref.resolveSibling("${ref.baseName}.dict"), checkIfExists: true)
    ch_assembly_ref = channel.value([ref, [ref_fai, ref_dict]])
    ch_ref_dict     = channel.value(ref_dict)

    ch_known_variants = channel.value([
        file(params.known_variants,          checkIfExists: true),
        file("${params.known_variants}.tbi", checkIfExists: true),
    ])
    ch_chunks_dir = channel.value(file(params.chunks_dir, checkIfExists: true))
    ch_glimpse_ref = channel.value([
        file(params.glimpse_ref_panel,          checkIfExists: true),
        file("${params.glimpse_ref_panel}.csi", checkIfExists: true),
    ])

    ch_highpass_vcf = channel.value([
        file(params.highpass_snp_vcf,          checkIfExists: true),
        file("${params.highpass_snp_vcf}.tbi", checkIfExists: true),
    ])

    // ── Shared VCF prep ──

    highpass_pass = FILTER_PASS_VCF(ch_highpass_vcf.map { vcf, tbi -> ["highpass", vcf, tbi] })

    ch_lowpass_vcfs_for_concat = channel.fromPath("${params.lowpass_qc_vcf_dir}/*.vcf.gz", checkIfExists: true)
        .collect().map { chr_vcf_files -> ["lowpass", chr_vcf_files] }
    ch_axiom_vcfs_for_concat = channel.fromPath("${params.axiom_lowpass_vcf_dir}/*.vcf.gz", checkIfExists: true)
        .collect().map { chr_vcf_files -> ["axiom", chr_vcf_files] }

    concat_out = CONCAT_CHR_VCFS(ch_lowpass_vcfs_for_concat.mix(ch_axiom_vcfs_for_concat))
    lowpass_allchr = concat_out.vcf.filter { label, _vcf, _idx -> label == "lowpass" }
    axiom_allchr   = concat_out.vcf.filter { label, _vcf, _idx -> label == "axiom" }

    ch_lowpass_11dogs_extract_in = lowpass_allchr
        .map { _label, vcf, idx -> ["lowpass_11dogs", vcf, idx] }
        .combine(channel.value(file(params.samples_highpass_11, checkIfExists: true)))
    ch_highpass_10dogs_extract_in = highpass_pass.vcf
        .map { _label, vcf, idx -> ["highpass_10dogs_down", vcf, idx] }
        .combine(channel.value(file(params.samples_downsampling_10, checkIfExists: true)))

    extract_out = EXTRACT_SAMPLES_VCF(ch_lowpass_11dogs_extract_in.mix(ch_highpass_10dogs_extract_in))
    lowpass_11       = extract_out.vcf.filter { label, _vcf, _idx -> label == "lowpass_11dogs" }
    truth_10dogs_raw = extract_out.vcf.filter { label, _vcf, _idx -> label == "highpass_10dogs_down" }

    ch_highpass_for_maf      = highpass_pass.vcf.map { _label, vcf, idx -> ["highpass", vcf, idx] }
    ch_lowpass_11_for_maf    = lowpass_11.map        { _label, vcf, idx -> ["lowpass_11", vcf, idx] }
    ch_lowpass_axiom_for_maf = lowpass_allchr.map    { _label, vcf, idx -> ["lowpass_axiom", vcf, idx] }
    ch_axiom_for_maf         = axiom_allchr.map      { _label, vcf, idx -> ["axiom", vcf, idx] }

    maf_out = MAF_FILTER_VCF(
        ch_highpass_for_maf.mix(ch_lowpass_11_for_maf, ch_lowpass_axiom_for_maf, ch_axiom_for_maf)
    )
    highpass_maf      = maf_out.vcf.filter { label, _vcf, _idx -> label == "highpass" }
    lowpass_11_maf    = maf_out.vcf.filter { label, _vcf, _idx -> label == "lowpass_11" }
    lowpass_axiom_maf = maf_out.vcf.filter { label, _vcf, _idx -> label == "lowpass_axiom" }
    axiom_maf         = maf_out.vcf.filter { label, _vcf, _idx -> label == "axiom" }

    ch_truth_10dogs = truth_10dogs_raw.map { _label, vcf, idx -> [vcf, idx] }

    // ── Downsampling study ──

    def FRACTIONS = ['0.1x', '0.2x', '0.3x', '0.4x', '0.5x', '0.6x', '0.7x', '0.8x', '0.9x', '1x']

    ch_samples_down = channel.fromPath(params.samples_downsampling_10, checkIfExists: true)
        .splitText()
        .map { line -> line.trim() }
        .filter { sample_id -> sample_id }

    ch_p_at_1x = channel.fromList(
        samplesheetToList(params.downsample_p_at_1x, "${projectDir}/assets/schemas/schema_downsample_p_at_1x.json")
    )

    ch_highpass_bam = ch_samples_down.map { sample_id ->
        [
            sample_id,
            file("${params.highpass_bam_dir}/${sample_id}.sorted.merged.MarkDups.BQSR.bam",     checkIfExists: true),
            file("${params.highpass_bam_dir}/${sample_id}.sorted.merged.MarkDups.BQSR.bam.bai", checkIfExists: true),
        ]
    }

    ch_downsample_in = ch_highpass_bam
        .combine(ch_p_at_1x, by: 0)
        .combine(channel.fromList(FRACTIONS))
        .map { sample_id, bam, bai, p_at_1x, fraction -> [sample_id, fraction, bam, bai, p_at_1x] }
    downsampled = DOWNSAMPLE_BAM(ch_downsample_in)

    ch_remap_in = downsampled.bam
        .map { sample_id, fraction, bam, bai -> [fraction, sample_id, [bam], [bai]] }
    remapped = MERGE_BQSR_GVCF_DOWN(ch_remap_in, ch_assembly_ref, ch_known_variants, ch_chunks_dir)

    ch_bam_bundles_by_fraction = remapped.bam
        .map { fraction, _sample_id, bam, bai -> [fraction, bam, bai] }
        .groupTuple()

    // Chunk/site definitions depend only on the reference panel, not on any sample or fraction.
    ch_chromosomes = channel.fromList((1..38).collect { num -> "chr${num}" })
    chunks_out = GLIMPSE_DEFINE_CHUNKS(ch_chromosomes, ch_glimpse_ref)
    sites_out  = GLIMPSE_VARIABLE_SITES(ch_chromosomes, ch_glimpse_ref)

    imputed = IMPUTE_DOWNSAMPLED_FRACTIONS(
        ch_bam_bundles_by_fraction, chunks_out.chunks, sites_out.sites, ch_assembly_ref, ch_glimpse_ref
    )
    // imputed.vcf: [label = "downsampled_<fraction>", vcf, csi]

    // ── Concordance comparisons ──

    ch_study1_comparison = lowpass_11_maf.map { _label, vcf, idx -> ["highpass_lowpass", vcf, idx] }
        .combine(highpass_maf.map { _label, vcf, idx -> [vcf, idx] })
    ch_study2_comparison = lowpass_axiom_maf.map { _label, vcf, idx -> ["axiom_lowpass", vcf, idx] }
        .combine(axiom_maf.map { _label, vcf, idx -> [vcf, idx] })
    ch_downsampled_comparisons = imputed.vcf.combine(ch_truth_10dogs)

    comparisons = ch_study1_comparison.mix(ch_study2_comparison, ch_downsampled_comparisons)

    pairs_study1 = channel.fromList(
        samplesheetToList(params.concordance_highpass_lowpass_pairs, "${projectDir}/assets/schemas/schema_concordance_pairs.json")
    ).map { pair_id, call_sample, truth_sample -> ["highpass_lowpass", pair_id, call_sample, truth_sample] }
    pairs_study2 = channel.fromList(
        samplesheetToList(params.concordance_axiom_lowpass_pairs, "${projectDir}/assets/schemas/schema_concordance_pairs.json")
    ).map { pair_id, call_sample, truth_sample -> ["axiom_lowpass", pair_id, call_sample, truth_sample] }
    // call_sample == truth_sample: same dog, compared against its own un-downsampled HighPass truth genotype.
    pairs_downsampling = ch_samples_down
        .combine(channel.fromList(FRACTIONS))
        .map { sample_id, fraction -> ["downsampled_${fraction}", "${sample_id}_${fraction}", sample_id, sample_id] }

    pairs = pairs_study1.mix(pairs_study2, pairs_downsampling)

    concordance_out = CONCORDANCE_ANALYSIS(comparisons, pairs, ch_ref_dict)

    ch_all_versions = channel.topic("versions")
        .unique()
        .map { process, tool, version ->
            [ process.tokenize(":").last(), "  ${tool}: ${version}" ]
        }
        .groupTuple(by: 0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }
        .collectFile(name: 'versions.yml', sort: true, newLine: true)

    publish:
    summary_metrics     = concordance_out.summary_metrics
    contingency_metrics = concordance_out.contingency_metrics
    versions            = ch_all_versions

}

// Output publishing mode is symlink by default. This saves space over copying.
// Filenames (e.g. "downsampled_0.1x_SRR..._concordance.vcf.genotype_concordance_summary_metrics")
// already identify which of the three studies each file belongs to.
output {
    summary_metrics {
        path 'concordance/summary_metrics'
    }
    contingency_metrics {
        path 'concordance/contingency_metrics'
    }
    versions {
        path 'pipeline_info'
    }
}
