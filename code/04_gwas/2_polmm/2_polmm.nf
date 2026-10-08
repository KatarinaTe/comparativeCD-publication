// ============================================================================
// main.nf — 04. GWAS — 2. POLMM workflow
//
// Per-item POLMM null model fit + full-genome marker test, for all 14 dogCD survey items, then
// plink clumping + sumstats export for all 14.
//
// Source: 04_gwas/2_polmm/{POLMMgrab_PER_ITEM.sh, clump_plink_items.sh,
//          export_gwas_sumstats_CCD_POLMM.R}
// ============================================================================

include { validateParameters } from 'plugin/nf-schema'

include { RUN_POLMM_ITEM      } from './process_definitions.nf'
include { CLUMP_ITEM          } from './process_definitions.nf'
include { EXPORT_SUMSTATS_POLMM } from './process_definitions.nf'

// The 14 dogCD survey items — matches 3_Merging_filtering's own getItemIds().
def getItemIds() {
    ['7', '93', '95', '145', '146', '147', '148', '149', '150', '151', '152', '153', '154', '155']
}

workflow {

    main:
    validateParameters()

    ch_gene_range_file = channel.value(file(params.gene_range_file, checkIfExists: true))

    ch_item_ids = getItemIds()

    ch_item_grab = channel.fromList(ch_item_ids)
        .map { id -> [id, file("${params.phenotypes_dir}/grabCCDitem${id}_DA_MERGED_GENCOVE_AXIOM_QC6.txt", checkIfExists: true)] }

    ch_item_qc6 = channel.fromList(ch_item_ids)
        .map { id ->
            def prefix = "CCDitem${id}_DA_MERGED_GENCOVE_AXIOM_QC6"
            [id,
             file("${params.qc6_plink_dir}/${prefix}.bed", checkIfExists: true),
             file("${params.qc6_plink_dir}/${prefix}.bim", checkIfExists: true),
             file("${params.qc6_plink_dir}/${prefix}.fam", checkIfExists: true)]
        }

    ch_item_eigeninput = channel.fromList(ch_item_ids)
        .map { id ->
            def prefix = "CCDitem${id}_EigenInput"
            [id,
             file("${params.qc6_plink_dir}/${prefix}.bed", checkIfExists: true),
             file("${params.qc6_plink_dir}/${prefix}.bim", checkIfExists: true),
             file("${params.qc6_plink_dir}/${prefix}.fam", checkIfExists: true)]
        }

    ch_polmm_input = ch_item_grab.join(ch_item_qc6).join(ch_item_eigeninput)

    polmm_out = RUN_POLMM_ITEM(ch_polmm_input)

    // ── Clumping + sumstats export — all 14 items ───────────────────────────

    ch_clump_input = ch_item_qc6.join(polmm_out.modi_sumstats)
    clump_out = CLUMP_ITEM(ch_clump_input, ch_gene_range_file)

    export_out = EXPORT_SUMSTATS_POLMM(polmm_out.modi_sumstats)

    publish:
    raw_marker        = polmm_out.raw_marker
    modi_sumstats     = polmm_out.modi_sumstats.map { _id, file -> file }
    clumped           = clump_out.clumped
    clumped_ranges    = clump_out.clumped_ranges
    clump_log         = clump_out.log
    sumstats          = export_out.sumstats

}

// Output publishing mode is symlink by default. This saves space over copying.
output {
    raw_marker {
        path 'gwas'
    }
    modi_sumstats {
        path 'gwas'
    }
    clumped {
        path 'clumping'
    }
    clumped_ranges {
        path 'clumping'
    }
    clump_log {
        path 'clumping'
    }
    sumstats {
        path 'sumstats'
    }
}
