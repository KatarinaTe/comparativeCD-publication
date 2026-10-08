// ============================================================================
// main.nf — 04. GWAS — 1. MLMA-LOCO workflow
//
// Per-phenotype mlma-loco GWAS, for the 3 CCD factors and SIZE (the control phenotype) — then
// plink clumping + sumstats export for all 4 phenotypes.
//
// SIZE is kept in its own published subdirectories from gene-finding (clumping) onward — the
// clumping step is where genes first get attached to a signal (--clump-range), and the sumstats
// files feed downstream finemapping/region-storage. Raw mlma-loco GWAS output (before
// gene-finding starts) stays in one shared directory across all 4 phenotypes.
//
// Source: 04_gwas/1_mlma-loco/{MLMA_PER_FACTOR.sh, clump_plink_factors.sh,
//          export_gwas_sumstats_CCD_MLMA.R}
// ============================================================================

include { validateParameters } from 'plugin/nf-schema'

include { RUN_MLMA_LOCO       } from './process_definitions.nf'
include { CLUMP_FACTOR        } from './process_definitions.nf'
include { EXPORT_SUMSTATS     } from './process_definitions.nf'

// F1/F2/F3 + SIZE all go through mlma-loco, clumping, and sumstats export uniformly — SIZE's
// outputs are kept in separate published subdirectories from clumping onward (see this file's
// header), not excluded from the processes themselves.
def getAllPhenoIds() { ['F1', 'F2', 'F3', 'SIZE'] }

// SIZE's QC6/phenotype files have no "CCD" prefix, unlike the factors (see
// 3_Merging_filtering's FILTER_SIZE_QC6/BUILD_SIZE_COVARIATES header comments for why).
def bfilePrefix(pheno_id) {
    pheno_id == 'SIZE' ? 'SIZE_DA_MERGED_GENCOVE_AXIOM_QC6' : "CCD${pheno_id}_DA_MERGED_GENCOVE_AXIOM_QC6"
}

workflow {

    main:
    validateParameters()

    ch_gene_range_file = channel.value(file(params.gene_range_file, checkIfExists: true))

    ch_pheno_plink = channel.fromList(getAllPhenoIds())
        .map { id ->
            def prefix = bfilePrefix(id)
            [id,
             file("${params.qc6_plink_dir}/${prefix}.bed", checkIfExists: true),
             file("${params.qc6_plink_dir}/${prefix}.bim", checkIfExists: true),
             file("${params.qc6_plink_dir}/${prefix}.fam", checkIfExists: true)]
        }

    ch_pheno_covariates = channel.fromList(getAllPhenoIds())
        .map { id ->
            def prefix = bfilePrefix(id)
            [id,
             file("${params.phenotypes_dir}/${prefix}.phen", checkIfExists: true),
             file("${params.phenotypes_dir}/age${prefix}.qcovar", checkIfExists: true),
             file("${params.phenotypes_dir}/sex${prefix}.covar", checkIfExists: true)]
        }

    // ── mlma-loco — all 4 phenotypes ────────────────────────────────────────

    ch_mlma_input = ch_pheno_plink.join(ch_pheno_covariates)
    mlma_out = RUN_MLMA_LOCO(ch_mlma_input)

    // ── Clumping + sumstats export — all 4 phenotypes ───────────────────────

    ch_clump_input = ch_pheno_plink.join(mlma_out.mlma)
    clump_out = CLUMP_FACTOR(ch_clump_input, ch_gene_range_file)

    ch_export_input = mlma_out.mlma
        .join(ch_pheno_covariates.map { id, phen, qcovar, covar -> [id, phen] })

    export_out = EXPORT_SUMSTATS(ch_export_input)

    // SIZE's clumped regions/sumstats are split into their own subdirectories from here on —
    // see this file's header. "SIZE_"-prefixed filenames vs. "CCD"-prefixed ones (set in
    // process_definitions.nf's CLUMP_FACTOR/EXPORT_SUMSTATS) tell them apart.
    ch_clumped_size    = clump_out.clumped.filter        { it.name.startsWith('SIZE_') }
    ch_clumped_factor  = clump_out.clumped.filter        { !it.name.startsWith('SIZE_') }
    ch_ranges_size     = clump_out.clumped_ranges.filter { it.name.startsWith('SIZE_') }
    ch_ranges_factor   = clump_out.clumped_ranges.filter { !it.name.startsWith('SIZE_') }
    ch_clumplog_size   = clump_out.log.filter            { it.name.startsWith('SIZE_') }
    ch_clumplog_factor = clump_out.log.filter            { !it.name.startsWith('SIZE_') }
    ch_sumstats_size   = export_out.sumstats.filter      { it.name.startsWith('SIZE_') }
    ch_sumstats_factor = export_out.sumstats.filter      { !it.name.startsWith('SIZE_') }

    publish:
    mlma_loco              = mlma_out.mlma.map { _id, mlma -> mlma }
    mlma_loco_log          = mlma_out.log
    clumped                = ch_clumped_factor
    clumped_size           = ch_clumped_size
    clumped_ranges         = ch_ranges_factor
    clumped_ranges_size    = ch_ranges_size
    clump_log              = ch_clumplog_factor
    clump_log_size         = ch_clumplog_size
    sumstats               = ch_sumstats_factor
    sumstats_size          = ch_sumstats_size

}

// Output publishing mode is symlink by default. This saves space over copying.
output {
    mlma_loco {
        path 'gwas'
    }
    mlma_loco_log {
        path 'gwas'
    }
    clumped {
        path 'clumping'
    }
    clumped_size {
        path 'clumping/size'
    }
    clumped_ranges {
        path 'clumping'
    }
    clumped_ranges_size {
        path 'clumping/size'
    }
    clump_log {
        path 'clumping'
    }
    clump_log_size {
        path 'clumping/size'
    }
    sumstats {
        path 'sumstats'
    }
    sumstats_size {
        path 'sumstats/size'
    }
}
