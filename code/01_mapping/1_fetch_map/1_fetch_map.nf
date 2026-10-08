// ============================================================================
// main.nf — 01_mapping / 1_fetch_map
//
// Per-run, per-sample-independent stage: FASTQ -> mapped/merged/BQSR'd BAM +
// gVCF + per-sample stats. Batchable — safe to invoke repeatedly over disjoint
// subsets of samples (see analyses/00_fetch-raw-data/run_batched.sh), unlike
// 2_cohort_genotyping which must only ever run once over the complete cohort.
//
// Track-agnostic by design: HighPass and LowPass samples go through exactly
// the same processes here. Which track a given invocation processes is
// determined by which samplesheet it's pointed at (params.samplesheet), not
// by separate process aliases — a natural consequence of a batch always being
// single-track. See archive/conversion_notes/01_mapping.md for the full
// rationale behind the two-pipeline split.
// ============================================================================

include { samplesheetToList; validateParameters } from 'plugin/nf-schema'

include { BWA_MEM2_INDEX                     } from './process_definitions.nf'
include { DOG10K_MAPPING                     } from './process_definitions.nf'
include { DOG10K_MERGE_MARKDUPS_BQSR_GVCF    } from './process_definitions.nf'
include { RUN_STATS                          } from './process_definitions.nf'

workflow {

    main:
    // ── Reference value channels ──────────────────────────────────────────────
    validateParameters()

    def ref      = file(params.assembly_ref,                        checkIfExists: true)
    def ref_fai  = file(ref.resolveSibling("${ref.name}.fai"),      checkIfExists: true)
    def ref_dict = file(ref.resolveSibling("${ref.baseName}.dict"), checkIfExists: true)
    index_out       = BWA_MEM2_INDEX(channel.value(ref))
    ch_assembly_ref = index_out.index
        .combine(channel.value(ref_fai))
        .combine(channel.value(ref_dict))
        .map { fasta, idx_files, fai, dict -> [fasta, [fai, dict] + idx_files] }
        .collect()
    ch_known_variants = channel.value([
        file(params.known_variants,          checkIfExists: true),
        file("${params.known_variants}.tbi", checkIfExists: true),
    ])
    ch_depth_sites = channel.value([
        file(params.depth_sites,             checkIfExists: true),
        file("${params.depth_sites}.tbi",    checkIfExists: true),
    ])
    ch_chunks_dir = channel.value(file(params.chunks_dir, checkIfExists: true))

    // ── Batch samplesheet: one track, one batch, per invocation ────────────────
    ch_reads = channel.fromList(samplesheetToList(
        params.samplesheet,
        "${projectDir}/assets/schemas/schema_input.json"
    ))

    mapping_out   = DOG10K_MAPPING(ch_reads, ch_assembly_ref)
    ch_per_sample = mapping_out.bam.groupTuple()

    merge_out = DOG10K_MERGE_MARKDUPS_BQSR_GVCF(
        ch_per_sample,
        ch_assembly_ref,
        ch_known_variants,
        ch_chunks_dir
    )

    // Per-sample stats — only the depth file is consumed downstream (by
    // 2_cohort_genotyping's EXTRACT_DEPTH, LowPass only), but it is produced
    // and materialized for both tracks, same as RUN_STATS always did.
    ch_for_stats = merge_out.bam
        .join(merge_out.gvcf)
        .map { sample_id, bam, bai, gvcf, gvcf_tbi ->
            [sample_id, bam, bai, gvcf, gvcf_tbi]
        }
    stats_out = RUN_STATS(ch_for_stats, ch_assembly_ref, ch_depth_sites)

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
    bam      = merge_out.bam
    gvcf     = merge_out.gvcf
    depth    = stats_out.depth
    versions = ch_all_versions

}

// Output publishing mode is symlink by default, same as the original combined
// pipeline — these paths still live inside this run's work/ directory.
// analyses/00_fetch-raw-data/run_batched.sh is what turns them into permanent,
// disk-safe files: it dereferences these symlinks into a permanent archive
// (cp -aL / rsync -L) before deleting this batch's work directory.
output {
    bam {
        path 'bam'
    }
    gvcf {
        path 'gvcf'
    }
    depth {
        path 'depth'
    }
    versions {
        path 'pipeline_info'
    }
}
