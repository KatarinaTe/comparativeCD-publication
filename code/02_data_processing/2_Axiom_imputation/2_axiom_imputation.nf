// ============================================================================
// main.nf — 2. Axiom imputation workflow
//
// Axiom array genotypes (canFam3) -> liftover to canFam4 -> per-chromosome
// imputation against the Dog10K phased panel -> merged PLINK dataset, plus the per-chromosome
// QC'd VCF (needed in VCF form, not just PLINK, by 03_qc's Axiom-vs-LowPass concordance study).
//
// Source scripts:
//   02_data_processing/2_Axiom_imputation/01_axiom_liftover.sh
//   02_data_processing/2_Axiom_imputation/02_axiom_impute.sh
//   02_data_processing/2_Axiom_imputation/03_axiomimp_to_plink.sh
// ============================================================================

include { validateParameters } from 'plugin/nf-schema'

include { LIFTOVER_TO_CANFAM4             } from './process_definitions.nf'
include { SPLIT_REF_PANEL_BY_CHR          } from './process_definitions.nf'
include { RECODE_AND_RENAME_CHR_VCF       } from './process_definitions.nf'
include { CONFORM_GT                      } from './process_definitions.nf'
include { BEAGLE_IMPUTE                   } from './process_definitions.nf'
include { QC_AND_CONVERT_CHR              } from './process_definitions.nf'
include { MERGE_CHR_PLINK_AXIOM           } from './process_definitions.nf'

workflow {

    main:
    validateParameters()

    // ── Reference value channels ────────────────────────────────────────────
    def ref     = file(params.assembly_ref,                   checkIfExists: true)
    def ref_fai = file(ref.resolveSibling("${ref.name}.fai"), checkIfExists: true)
    ch_assembly_ref = channel.value([ref, ref_fai])

    ch_glimpse_ref = channel.value([
        file(params.glimpse_ref_panel,              checkIfExists: true),
        file("${params.glimpse_ref_panel}.csi",     checkIfExists: true),
    ])

    ch_ind_to_keep    = channel.value(file(params.ind_to_keep,     checkIfExists: true))
    ch_chain          = channel.value(file(params.liftover_chain,  checkIfExists: true))
    ch_chr_rename_map = channel.value(file(params.chr_rename_map,  checkIfExists: true))
    ch_conform_gt_jar = channel.value(file(params.conform_gt_jar,  checkIfExists: true))
    ch_beagle_jar     = channel.value(file(params.beagle_jar,      checkIfExists: true))

    ch_axiom = channel.value([
        file("${params.axiom_bfile_prefix}.bed", checkIfExists: true),
        file("${params.axiom_bfile_prefix}.bim", checkIfExists: true),
        file("${params.axiom_bfile_prefix}.fam", checkIfExists: true),
    ])

    // Per-chromosome [chr, start, stop] — drives ref-panel splitting and imputation
    ch_chr_start_stop = channel.fromPath(params.chr_start_stop, checkIfExists: true)
        .splitText()
        .map { line ->
            def fields = line.trim().split(/\s+/)
            [fields[0], fields[1], fields[2]]
        }

    // ── Liftover chain: canFam3 Axiom plinkset -> canFam4 (DA_AFFY) ─────────

    liftover_out = LIFTOVER_TO_CANFAM4(ch_axiom, ch_ind_to_keep, ch_chain)

    // Broadcast the single cohort-level DA_AFFY plinkset across all 38 per-chr emissions
    ch_da_affy = liftover_out.plink.collect()

    // ── Per-chromosome imputation against the Dog10K reference panel ───────

    ref_chr_out = SPLIT_REF_PANEL_BY_CHR(ch_chr_start_stop, ch_glimpse_ref)

    recode_rename_out = RECODE_AND_RENAME_CHR_VCF(ch_chr_start_stop, ch_da_affy, ch_assembly_ref, ch_chr_rename_map)

    ch_ref_by_chr = ref_chr_out.ref_chr
        .map { chr, _start, _stop, ref_vcf, ref_csi -> [chr, ref_vcf, ref_csi] }

    ch_for_conform = recode_rename_out.vcf
        .join(ch_ref_by_chr)
    // [chr, start, stop, gt_vcf, gt_csi, ref_vcf, ref_csi]

    conform_out = CONFORM_GT(ch_for_conform, ch_conform_gt_jar)

    ch_for_beagle = conform_out.vcf
        .join(ch_ref_by_chr)
        .map { chr, start, stop, conf_vcf, conf_csi, ref_vcf, ref_csi ->
            def genetic_map = file("${params.genetic_maps_dir}/canFam4.cM.${chr}.map", checkIfExists: true)
            [chr, start, stop, conf_vcf, conf_csi, ref_vcf, ref_csi, genetic_map]
        }

    beagle_out = BEAGLE_IMPUTE(ch_for_beagle, ch_beagle_jar)

    // ── QC filter, ID renaming, diagnostic counts, and PLINK conversion, then merge ─

    plink_out = QC_AND_CONVERT_CHR(beagle_out.vcf)

    ch_all_plink = plink_out.plink
        .flatMap { _chr, bed, bim, fam -> [bed, bim, fam] }
        .collect()
    merge_out = MERGE_CHR_PLINK_AXIOM(ch_all_plink)

    // Per-chr QC'd VCF, published alongside the merged PLINK conversion — needed by 03_qc's
    // Axiom-vs-LowPass concordance study (see QC_AND_CONVERT_CHR's header comment).
    ch_per_chr_vcf = plink_out.vcf.map { _chr, vcf, csi -> [vcf, csi] }.flatten()

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
    axiom_plink_bed = merge_out.bed
    axiom_plink_bim = merge_out.bim
    axiom_plink_fam = merge_out.fam
    axiom_plink_log = merge_out.log
    axiom_per_chr_vcf = ch_per_chr_vcf
    versions        = ch_all_versions

}

// Output publishing mode is symlink by default. This saves space over copying.
output {
    axiom_plink_bed {
        path 'plink'
    }
    axiom_plink_bim {
        path 'plink'
    }
    axiom_plink_fam {
        path 'plink'
    }
    axiom_plink_log {
        path 'plink'
    }
    axiom_per_chr_vcf {
        path 'vcf'
    }
    versions {
        path 'pipeline_info'
    }
}
