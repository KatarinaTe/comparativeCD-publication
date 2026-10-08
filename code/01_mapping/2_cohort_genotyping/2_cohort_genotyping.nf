// ============================================================================
// main.nf — 01_mapping / 2_cohort_genotyping
//
// Cohort-wide stage: joint SNP calling (HighPass) + GLIMPSE imputation
// (LowPass). MUST be invoked exactly ONCE, over the complete accumulated
// sample set, after every batch of 1_fetch_map has finished and been
// materialized for both tracks (see analyses/00_fetch-raw-data/run_batched.sh).
//
// This is a hard correctness constraint, not a performance one:
//   - GENOTYPE_CALLING joint-calls all HighPass gVCFs together. Splitting
//     that across batches would change the actual variant calls.
//   - GLIMPSE_LIKELIHOODS computes genotype likelihoods from the bundle of
//     ALL LowPass BAMs at once. GLIMPSE's imputation accuracy depends on
//     genotype likelihoods computed jointly across the whole cohort — running
//     this over a subset of samples would silently produce different, less
//     accurate imputation, not just run faster.
//
// Inputs are two samplesheets (HighPass/LowPass) of already-produced,
// already-merged per-sample BAM/GVCF/depth files — the same tuple shape that
// 1_fetch_map's DOG10K_MERGE_MARKDUPS_BQSR_GVCF/RUN_STATS emit, now sourced
// from files on disk (materialized by the batch driver) instead of an
// in-memory channel from the same workflow run. See
// assets/schemas/schema_cohort_input.json.
// ============================================================================

include { samplesheetToList; validateParameters } from 'plugin/nf-schema'

include { BWA_MEM2_INDEX          } from './process_definitions.nf'
include { GENOTYPE_CALLING        } from './process_definitions.nf'
include { EXTRACT_DEPTH           } from './process_definitions.nf'
include { GLIMPSE_DEFINE_CHUNKS   } from './process_definitions.nf'
include { GLIMPSE_VARIABLE_SITES  } from './process_definitions.nf'
include { GLIMPSE_LIKELIHOODS     } from './process_definitions.nf'
include { GLIMPSE_IMPUTE_CHUNKS   } from './process_definitions.nf'
include { INDEX_BCFS              } from './process_definitions.nf'
include { GLIMPSE_LIGATE          } from './process_definitions.nf'
include { BCF_VCF                 } from './process_definitions.nf'
include { GLIMPSE_QC_FILTER       } from './process_definitions.nf'
include { VCF_TO_PLINK            } from './process_definitions.nf'
include { MERGE_CHR_PLINK         } from './process_definitions.nf'

workflow {

    main:
    // ── Reference value channels ──────────────────────────────────────────────
    validateParameters()

    // Rebuilt the same way 1_fetch_map builds it (BWA_MEM2_INDEX is not
    // actually needed by any process below — GATK/bcftools only use .fai/
    // .dict — but the process input signatures are unchanged from the
    // original pipeline, which bundled the bwa-mem2 index files into the same
    // ch_assembly_ref tuple everywhere). See archive/conversion_notes/01_mapping.md.
    def ref      = file(params.assembly_ref,                        checkIfExists: true)
    def ref_fai  = file(ref.resolveSibling("${ref.name}.fai"),      checkIfExists: true)
    def ref_dict = file(ref.resolveSibling("${ref.baseName}.dict"), checkIfExists: true)
    index_out       = BWA_MEM2_INDEX(channel.value(ref))
    ch_assembly_ref = index_out.index
        .combine(channel.value(ref_fai))
        .combine(channel.value(ref_dict))
        .map { fasta, idx_files, fai, dict -> [fasta, [fai, dict] + idx_files] }
        .collect()
    ch_chunks_dir = channel.value(file(params.chunks_dir, checkIfExists: true))

    // GLIMPSE reference panel (Dog10K phased BCF + index)
    ch_glimpse_ref = channel.value([
        file(params.glimpse_ref_panel,          checkIfExists: true),
        file("${params.glimpse_ref_panel}.csi", checkIfExists: true),
    ])

    // Dog autosomes — drives all per-chromosome GLIMPSE processes
    ch_chromosomes = channel.fromList((1..38).collect { num -> "chr${num}" })

    // ── Track 1: HighPass — joint genotype calling ────────────────────────────

    ch_hp_samples = channel.fromList(samplesheetToList(
        params.high_pass_samplesheet,
        "${projectDir}/assets/schemas/schema_cohort_input.json"
    ))

    ch_hp_gvcfs = ch_hp_samples
        .map { _sample, _bam, _bai, gvcf, _gvcf_tbi, _depth -> gvcf }
        .collect()
    ch_hp_tbis = ch_hp_samples
        .map { _sample, _bam, _bai, _gvcf, gvcf_tbi, _depth -> gvcf_tbi }
        .collect()

    genotype_out = GENOTYPE_CALLING(
        ch_hp_gvcfs,
        ch_hp_tbis,
        ch_assembly_ref,
        ch_chunks_dir
    )

    // ── Track 2: LowPass — cohort depth + GLIMPSE imputation ──────────────────

    ch_lp_samples = channel.fromList(samplesheetToList(
        params.low_pass_samplesheet,
        "${projectDir}/assets/schemas/schema_cohort_input.json"
    ))

    // Cohort depth aggregation
    depth_stats = EXTRACT_DEPTH(
        ch_lp_samples
            .map { _sample, _bam, _bai, _gvcf, _gvcf_tbi, depth -> depth }
            .collect()
    )

    // Steps 1-2 run in parallel — both driven by ch_chromosomes
    chunks_out = GLIMPSE_DEFINE_CHUNKS(ch_chromosomes, ch_glimpse_ref)
    sites_out  = GLIMPSE_VARIABLE_SITES(ch_chromosomes, ch_glimpse_ref)

    // Collect all LowPass BAMs into a single bundle for genotype likelihoods.
    // toList() preserves tuple structure; we separate bams and bais.
    // combine() pairs each of 38 chr emissions with the one BAM-bundle emission.
    ch_lp_bam_bundle = ch_lp_samples
        .map { _sample, bam, bai, _gvcf, _gvcf_tbi, _depth -> [bam, bai] }
        .toList()
        .map { pairs ->
            def bams = pairs.collect { pair -> pair[0] }
            def bais = pairs.collect { pair -> pair[1] }
            [bams, bais]
        }

    ch_for_gl = sites_out.sites.combine(ch_lp_bam_bundle)
    // Each emission: [chr, sites_vcf, sites_csi, sites_tsv, sites_tbi, [bams], [bais]]

    // Step 3: compute genotype likelihoods per chromosome — ALL LowPass BAMs at once
    gl_out = GLIMPSE_LIKELIHOODS(ch_for_gl, ch_assembly_ref)

    // Step 4: impute each chunk independently
    // Split chunk definition files into individual records, then join with GL VCFs
    ch_chunk_records = chunks_out.chunks
        .flatMap { _chr, chunks_file ->
            chunks_file.readLines().collect { line ->
                def fields = line.trim().split(/\s+/)
                // GLIMPSE chunk format: ID CHR INPUT_REGION OUTPUT_REGION
                [fields[1], fields[0], fields[2], fields[3]]
                // [chr, chunk_id, input_region, output_region]
            }
        }

    ch_for_impute = ch_chunk_records
        .combine(gl_out.gl_vcf, by: 0)
    // [chr, chunk_id, input_region, output_region, gl_vcf, gl_vcf_csi]

    // Step 5a: index imputed BCF files (bcftools only)
    impute_out  = GLIMPSE_IMPUTE_CHUNKS(ch_for_impute, ch_glimpse_ref)
    indexed_out = INDEX_BCFS(impute_out.imputed_bcf)

    // Step 5b: ligate chunks per chromosome
    // Group indexed chunks by chr — sorting happens inside the process script
    ch_for_ligate = indexed_out.imputed
        .map { chr, _chunk_id, bcf, csi -> [chr, bcf, csi] }
        .groupTuple()

    ligate_out = GLIMPSE_LIGATE(ch_for_ligate)

    // Step 5c: convert merged BCF → bgzipped VCF (bcftools only)
    ch_for_bcf_vcf = ligate_out.merged_bcf
        .groupTuple()

    bcf_vcf_out = BCF_VCF(ch_for_bcf_vcf)

    // Step 6: QC filter per chromosome
    // Published (in addition to feeding VCF_TO_PLINK below) so downstream stages
    // (e.g. 03_qc concordance studies) can consume the low-pass cohort as a VCF
    // without round-tripping through the published PLINK files.
    qc_out = GLIMPSE_QC_FILTER(bcf_vcf_out.merged)

    // Step 7: convert to PLINK per chromosome
    plink_out = VCF_TO_PLINK(qc_out.qc_vcf)

    // Step 8: merge all chromosomes into final PLINK dataset
    // Flatten [chr, bed, bim, fam] → [bed, bim, fam] then collect all files
    ch_all_plink = plink_out.plink
        .flatMap { _chr, bed, bim, fam -> [bed, bim, fam] }
        .collect()
    merge_out = MERGE_CHR_PLINK(ch_all_plink)

    ch_all_versions = channel.topic("versions")
        .unique()
        .map { process, tool, version ->
            [ process.tokenize(":").last(), "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }
        .collectFile(name: 'versions.yml', sort: true, newLine: true)

    publish:
    lowpass_plink_bed = merge_out.bed
    lowpass_plink_bim = merge_out.bim
    lowpass_plink_fam = merge_out.fam
    lowpass_plink_log = merge_out.log
    lowpass_depth     = depth_stats.depth_stats
    highpass_snp_vcf  = genotype_out.vcf
    highpass_snp_tbi  = genotype_out.vcf_tbi
    lowpass_qc_vcf    = qc_out.qc_vcf
    versions          = ch_all_versions

}

// Output publishing mode is symlink by default. This saves space over copying.
// NOTE: highpass_bam/lowpass_bam are NOT published here — unlike the original
// combined 01_mapping.nf, per-sample BAMs are produced and materialized by
// 1_fetch_map (analyses/00_fetch-raw-data/run_batched.sh), upstream of this
// pipeline. Everything published below is genuinely downstream of per-sample
// BAM/GVCF, matching the original publish block for those outputs exactly.
output {
    lowpass_plink_bed {
        path { 'lowpass/plink' }
    }
    lowpass_plink_bim {
        path { 'lowpass/plink' }
    }
    lowpass_plink_fam {
        path { 'lowpass/plink' }
    }
    lowpass_plink_log {
        path { 'lowpass/plink' }
    }
    lowpass_depth {
        path 'lowpass/depth'
    }
    highpass_snp_vcf {
        path 'highpass/snp_vcf'
    }
    highpass_snp_tbi {
        path 'highpass/snp_vcf'
    }
    lowpass_qc_vcf {
        path 'lowpass/vcf'
    }
    versions {
        path 'pipeline_info'
    }
}
