// ============================================================================
// main.nf — 04. GWAS — 3. Finemap (SuSiE, no functional priors) workflow
//
// Builds one shared genofile (split per chromosome) that all phenotypes finemap against, then
// harmonizes + munges each of the 18 GWAS'd phenotypes' sumstats (3 factors + 14 items, plus
// SIZE), then runs PolyFun's finemapper.py for every one of the 113 windows in
// assets/finemap_regions.csv (104 hand-picked) + assets/finemap_regions_size.csv (9, mechanically
// extracted from SIZE's own already-run results). SIZE goes through the identical processes but
// is kept in its own published subdirectories throughout (sumstats/size, finemap/size).
//
// Source: 04_gwas/3_finemap_susie/{1_finemap_create_genofiles.sh,
//          2a_run_harmonize_sumstats_CCD.sbatch, harmonize_sumstats_CCD_SIZE.py,
//          2b_munge_sumstats_parquet.sh, finemapNOFUNCT_multiple_CCDregions{1,2,3,4}.sh},
//          data/04_gwas/SIZE_260603_CSnonfunct.xlsx (SIZE's regions; see extract_size_regions.py)
// ============================================================================

include { validateParameters } from 'plugin/nf-schema'

include { BUILD_COMBINED_FAM    } from './process_definitions.nf'
include { FILTER_COMBINED_QC6   } from './process_definitions.nf'
include { SPLIT_GENOFILE_BY_CHR } from './process_definitions.nf'
include { HARMONIZE_SUMSTATS    } from './process_definitions.nf'
include { MUNGE_SUMSTATS        } from './process_definitions.nf'
include { FINEMAP_REGION        } from './process_definitions.nf'

// Dog autosomes 1-38 — matches the original's own `for chr in {1..38}` loop.
def getChrList() { (1..38).collect { it as String } }

// The 17 phenotypes harmonized/munged/finemapped in the original scripts (STUCK excluded
// throughout this stage — see process_definitions.nf's header). SIZE is wired in separately
// below: data/04_gwas/SIZE_260603_CSnonfunct.xlsx is evidence SIZE was finemapped in the real
// analysis, with its own region windows recoverable from that file (see extract_size_regions.py).
def getFactorIds() { ['F1', 'F2', 'F3'] }
def getItemIds()   { ['7', '93', '95', '145', '146', '147', '148', '149', '150', '151', '152', '153', '154', '155'] }

workflow {

    main:
    validateParameters()

    // ── Shared genofile creation ─────────────────────────────────────────────

    // QC5 lives in the same published plink/ directory as QC6 (3_Merging_filtering publishes
    // both there) — reusing the same qc6_plink_dir param 1_mlma-loco/2_polmm already use, not a
    // separate qc5_plink_dir.
    ch_qc5_fam    = channel.value(file("${params.qc6_plink_dir}/DA_MERGED_GENCOVE_AXIOM_QC5.fam", checkIfExists: true))
    ch_qc5_bed    = channel.value(file("${params.qc6_plink_dir}/DA_MERGED_GENCOVE_AXIOM_QC5.bed", checkIfExists: true))
    ch_qc5_bim    = channel.value(file("${params.qc6_plink_dir}/DA_MERGED_GENCOVE_AXIOM_QC5.bim", checkIfExists: true))
    ch_size_fam   = channel.value(file(params.size_fam, checkIfExists: true))
    ch_stuck_fam  = channel.value(file(params.stuck_fam, checkIfExists: true))
    ch_allfam_fam = channel.value(file(params.allfam_fam, checkIfExists: true))

    combined_fam_out = BUILD_COMBINED_FAM(ch_qc5_fam, ch_size_fam, ch_stuck_fam, ch_allfam_fam)

    qc6_out = FILTER_COMBINED_QC6(ch_qc5_bed.combine(ch_qc5_bim), combined_fam_out.fam)

    ch_chr_input = channel.fromList(getChrList())
        .combine(qc6_out.plink)
        .map { chr, bed, bim, fam -> [chr, bed, bim, fam] }

    chr_geno_out = SPLIT_GENOFILE_BY_CHR(ch_chr_input)

    // ── Harmonize + munge — 18 phenotypes (3 factors + 14 items + SIZE) ─────

    ch_factor_sumstats = channel.fromList(getFactorIds())
        .map { id -> ["CCD${id}", file("${params.mlma_sumstats_dir}/CCD${id}_sumstats.txt", checkIfExists: true)] }

    ch_item_sumstats = channel.fromList(getItemIds())
        .map { id -> ["ITEM${id}", file("${params.polmm_sumstats_dir}/ITEM${id}_sumstats.txt", checkIfExists: true)] }

    // SIZE's sumstats live in their own subdirectory (1_mlma-loco's sumstats/size/), not the flat
    // sumstats/ dir the factors use — see 1_mlma-loco's own "keep SIZE separated" publish split.
    ch_size_sumstats = channel.of('SIZE')
        .map { id -> [id, file("${params.size_sumstats_dir}/SIZE_sumstats.txt", checkIfExists: true)] }

    ch_pheno_sumstats = ch_factor_sumstats.mix(ch_item_sumstats).mix(ch_size_sumstats)

    // .collect() converts this single-emission process-output channel into a broadcastable value
    // channel — without it, Nextflow zips it against the 17-element sumstats channel instead of
    // reusing it for every element, so HARMONIZE_SUMSTATS would only ever fire once. Preferred
    // over .first() here: .first() silently keeps only the first emission and drops the rest if
    // this channel ever unexpectedly emits more than one item, masking a real bug; .collect()
    // gathers everything, so an unexpected extra emission stays visible instead of being dropped.
    ch_combined_bim = qc6_out.plink.map { bed, bim, fam -> bim }.collect()

    harmonize_out = HARMONIZE_SUMSTATS(ch_pheno_sumstats, ch_combined_bim)
    munge_out     = MUNGE_SUMSTATS(harmonize_out.harmonized)

    // ── Finemap — 104 hand-picked + 9 SIZE windows ───────────────────────────

    ch_regions_main = channel.fromPath("${projectDir}/assets/finemap_regions.csv")
        .splitCsv(header: true)
    ch_regions_size = channel.fromPath("${projectDir}/assets/finemap_regions_size.csv")
        .splitCsv(header: true)
    ch_regions = ch_regions_main.mix(ch_regions_size)

    ch_regions_with_geno = ch_regions
        .map { row -> [row.chr, row] }
        .combine(chr_geno_out.plink, by: 0)
        .map { chr, row, bed, bim, fam -> [row.pheno_id, row, bed, bim, fam] }

    ch_finemap_input = ch_regions_with_geno
        .combine(munge_out.munged, by: 0)
        .map { pheno_id, row, bed, bim, fam, parquet -> [row, bed, bim, fam, parquet] }

    finemap_out = FINEMAP_REGION(ch_finemap_input)

    // SIZE is kept in its own published subdirectories from here on, same as 1_mlma-loco's own
    // "keep SIZE separated" split — it goes through the identical HARMONIZE_SUMSTATS/
    // MUNGE_SUMSTATS/FINEMAP_REGION processes, but its outputs never mix into the main
    // sumstats//finemap/ directories. All of SIZE's own output filenames start with "SIZE_"
    // (pheno_id "SIZE"), which no other of the 17 phenotypes' IDs are a prefix of.
    ch_harmonized_size    = harmonize_out.harmonized.filter   { id, _f -> id == 'SIZE' }.map { _id, f -> f }
    ch_harmonized_main    = harmonize_out.harmonized.filter   { id, _f -> id != 'SIZE' }.map { _id, f -> f }
    ch_swapped_size       = harmonize_out.swapped_snps.filter { it.name.startsWith('SIZE_') }
    ch_swapped_main       = harmonize_out.swapped_snps.filter { !it.name.startsWith('SIZE_') }
    ch_munged_size        = munge_out.munged.filter           { id, _f -> id == 'SIZE' }.map { _id, f -> f }
    ch_munged_main        = munge_out.munged.filter           { id, _f -> id != 'SIZE' }.map { _id, f -> f }
    ch_finemap_size       = finemap_out.result.filter         { it.name.startsWith('SIZE_') }
    ch_finemap_main       = finemap_out.result.filter         { !it.name.startsWith('SIZE_') }

    publish:
    combined_qc6           = qc6_out.plink.flatMap { bed, bim, fam -> [bed, bim, fam] }
    combined_qc6_log       = qc6_out.log
    chr_geno               = chr_geno_out.plink.flatMap { _chr, bed, bim, fam -> [bed, bim, fam] }
    harmonized_sumstats    = ch_harmonized_main
    harmonized_sumstats_size = ch_harmonized_size
    swapped_snps            = ch_swapped_main
    swapped_snps_size        = ch_swapped_size
    munged_sumstats          = ch_munged_main
    munged_sumstats_size     = ch_munged_size
    finemap_results          = ch_finemap_main
    finemap_results_size     = ch_finemap_size

}

// Output publishing mode is symlink by default. This saves space over copying.
output {
    combined_qc6 {
        path 'genofiles'
    }
    combined_qc6_log {
        path 'genofiles'
    }
    chr_geno {
        path 'genofiles/per_chr'
    }
    harmonized_sumstats {
        path 'sumstats'
    }
    harmonized_sumstats_size {
        path 'sumstats/size'
    }
    swapped_snps {
        path 'sumstats/diagnostics'
    }
    swapped_snps_size {
        path 'sumstats/size/diagnostics'
    }
    munged_sumstats {
        path 'sumstats'
    }
    munged_sumstats_size {
        path 'sumstats/size'
    }
    finemap_results {
        path 'finemap'
    }
    finemap_results_size {
        path 'finemap/size'
    }
}
