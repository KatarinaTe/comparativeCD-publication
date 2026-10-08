// ============================================================================
// main.nf — 3. Merging & filtering workflow
//
// Merges the LowPass GENCOVE (01_mapping) and Axiom (2_Axiom_imputation) genotype datasets,
// filters samples and variants, and produces per-phenotype (3 factors + 14 items) GWAS-ready
// PLINK datasets.
//
// Source: 02_data_processing/3_Merging_filtering/{01_gencove_axiom_merging_filtering.sh,
//          02_create_ALLFAM_check_stats.sh, 03_create_files_per_factor_mlma.R,
//          04_create_files_per_item_polmm.R}
//
// Implements Stage A (cohort QC + merge), Stage B (flip-scan convergence), Stage C (sample-level
// QC: duplicates, depth, KING relatedness, dogID relabeling), Stage D (ALLFAM stats, PCA, GRM),
// Stage E (per-factor phenotype files), Stage F (per-item phenotype files), and Stage G
// (cross-phenotype dog/SNP-count check).
//
// QC5 excludes 12 additional LowPass samples that 01_gencove_axiom_merging_filtering.sh never
// actually removes (see APPLY_EXCLUDE_LOWPASS12/FILTER_DOGID_RELABEL_FAM), shrinking it to 3316
// dogs. The frozen data6 checkpoint is untouched by this and still has 3328 rows, so it is
// restricted + reordered to QC5's own dogID set/order (FILTER_DATA6_TO_QC5) before
// BUILD_FACTOR_FAM/BUILD_ITEM_FAM use it.
// ============================================================================

include { validateParameters } from 'plugin/nf-schema'

include { QC_FILTER_GENCOVE        } from './process_definitions.nf'
include { QC_FILTER_AXIOM          } from './process_definitions.nf'
include { MERGE_DIAGNOSTIC         } from './process_definitions.nf'
include { MERGE_COHORTS            } from './process_definitions.nf'
include { QC_MERGED                } from './process_definitions.nf'
include { FLIP_SCAN_CONVERGENCE    } from './process_definitions.nf'
include { QC_FINAL_MERGE           } from './process_definitions.nf'
include { BUILD_EXCLUDE_LIST_ROUND1} from './process_definitions.nf'
include { APPLY_EXCLUDE_ROUND1     } from './process_definitions.nf'
include { PRUNE_FOR_KING           } from './process_definitions.nf'
include { PCA_ROUND1               } from './process_definitions.nf'
include { RUN_KING_ROUND1          } from './process_definitions.nf'
include { APPLY_EXCLUDE_ROUND2     } from './process_definitions.nf'
include { PCA_ROUND2               } from './process_definitions.nf'
include { RUN_KING_ROUND2          } from './process_definitions.nf'
include { APPLY_EXCLUDE_LOWPASS12  } from './process_definitions.nf'
include { FILTER_DOGID_RELABEL_FAM } from './process_definitions.nf'
include { RELABEL_AND_EXCLUDE_ROUND3 } from './process_definitions.nf'
include { BUILD_DATA1              } from './process_definitions.nf'
include { BUILD_SIZE_FAM           } from './process_definitions.nf'
include { FILTER_SIZE_QC6          } from './process_definitions.nf'
include { BUILD_SIZE_COVARIATES    } from './process_definitions.nf'
include { BUILD_ALLFAM_QC6         } from './process_definitions.nf'
include { PRUNE_ALLFAM             } from './process_definitions.nf'
include { PCA_ALLFAM               } from './process_definitions.nf'
include { BUILD_GRM                } from './process_definitions.nf'
include { DEPTH_HISTOGRAMS         } from './process_definitions.nf'
include { FILTER_DATA6_TO_QC5      } from './process_definitions.nf'
include { BUILD_FACTOR_FAM         } from './process_definitions.nf'
include { FILTER_FACTOR_QC6        } from './process_definitions.nf'
include { BUILD_FACTOR_COVARIATES  } from './process_definitions.nf'
include { BUILD_ITEM_FAM           } from './process_definitions.nf'
include { FILTER_ITEM_QC6          } from './process_definitions.nf'
include { PRUNE_ITEM_LD            } from './process_definitions.nf'
include { BUILD_ITEM_EIGENINPUT    } from './process_definitions.nf'
include { BUILD_ITEM_GRAB          } from './process_definitions.nf'
include { CHECK_PHENOTYPE_COUNTS   } from './process_definitions.nf'

// Fixed per-factor configuration — literal values from the original script (item columns feeding
// each factor's NA count, and each factor's NA-count exclusion threshold), not runtime input.
// See bin/build_factor_fam.R's header for the na_count bug fix this config feeds.
def getFactorConfigs() {
    [
        [factor_id: 'F1', item_cols: 'item7,item153,item154,item155', na_threshold: 3],
        [factor_id: 'F2', item_cols: 'item93,item95,item150',         na_threshold: 2],
        [factor_id: 'F3', item_cols: 'item145,item146,item147,item148,item149,item151,item152', na_threshold: 7],
    ]
}

// The 14 dogCD survey items, each getting its own per-item phenotype pipeline (Stage F).
def getItemIds() {
    ['7', '93', '95', '145', '146', '147', '148', '149', '150', '151', '152', '153', '154', '155']
}

workflow {

    main:
    validateParameters()

    ch_gencove = channel.value([
        file("${params.gencove_bfile_prefix}.bed", checkIfExists: true),
        file("${params.gencove_bfile_prefix}.bim", checkIfExists: true),
        file("${params.gencove_bfile_prefix}.fam", checkIfExists: true),
    ])
    ch_axiom = channel.value([
        file("${params.axiom_bfile_prefix}.bed", checkIfExists: true),
        file("${params.axiom_bfile_prefix}.bim", checkIfExists: true),
        file("${params.axiom_bfile_prefix}.fam", checkIfExists: true),
    ])

    ch_fam_qc3modi       = channel.value(file(params.fam_qc3modi,       checkIfExists: true))
    ch_depth_file        = channel.value(file(params.depth_file,        checkIfExists: true))
    ch_duplicates_file   = channel.value(file(params.duplicates_file,   checkIfExists: true))
    ch_exclude_round2    = channel.value(file(params.exclude_list_round2, checkIfExists: true))
    ch_dogid_relabel_fam = channel.value(file(params.dogid_relabel_fam, checkIfExists: true))
    ch_exclude_round3    = channel.value(file(params.exclude_list_round3, checkIfExists: true))
    ch_covariates        = channel.value(file(params.covariates,        checkIfExists: true))
    ch_fam_2584dogs      = channel.value(file(params.fam_2584dogs,      checkIfExists: true))
    ch_data6             = channel.value(file(params.data6,             checkIfExists: true))
    ch_exclude_lowpass12 = channel.value(file(params.exclude_list_lowpass12, checkIfExists: true))
    ch_size_pheno_file   = channel.value(file(params.size_pheno_file,    checkIfExists: true))

    // ── Stage A: cohort QC + merge ───────────────────────────────────────────

    gencove_qc_out = QC_FILTER_GENCOVE(ch_gencove)
    axiom_qc_out   = QC_FILTER_AXIOM(ch_axiom)

    MERGE_DIAGNOSTIC(gencove_qc_out.plink, axiom_qc_out.plink)
    merge_out = MERGE_COHORTS(gencove_qc_out.plink, axiom_qc_out.plink)
    qc_merged_out = QC_MERGED(merge_out.plink)

    // ── Stage B: flip-scan convergence ──────────────────────────────────────

    flipscan_out = FLIP_SCAN_CONVERGENCE(
        qc_merged_out.plink,
        gencove_qc_out.plink,
        axiom_qc_out.plink
    )

    // ── Stage C: sample-level QC — duplicates, depth, KING, dogID relabeling ─

    qc4_out = QC_FINAL_MERGE(flipscan_out.plink)

    exclude1_out = BUILD_EXCLUDE_LIST_ROUND1(
        ch_fam_qc3modi,
        ch_depth_file,
        ch_duplicates_file,
        qc4_out.plink
    )

    qc4b_out = APPLY_EXCLUDE_ROUND1(qc4_out.plink, exclude1_out.exclude_list)

    prune_out = PRUNE_FOR_KING(qc4b_out.plink)

    pca1_out = PCA_ROUND1(qc4b_out.plink, prune_out.prune_in)
    RUN_KING_ROUND1(pca1_out.plink)

    qc4c_out = APPLY_EXCLUDE_ROUND2(qc4b_out.plink, ch_exclude_round2)

    pca2_out = PCA_ROUND2(qc4c_out.plink, prune_out.prune_in)
    RUN_KING_ROUND2(pca2_out.plink)

    // Not in the original script — see APPLY_EXCLUDE_LOWPASS12's header comment.
    lowpass12_qc4c_out = APPLY_EXCLUDE_LOWPASS12(qc4c_out.plink, ch_exclude_lowpass12)
    lowpass12_fam_out  = FILTER_DOGID_RELABEL_FAM(ch_dogid_relabel_fam, ch_exclude_lowpass12)

    qc5_out = RELABEL_AND_EXCLUDE_ROUND3(
        lowpass12_qc4c_out.plink,
        lowpass12_fam_out.fam,
        ch_exclude_round3
    )

    data1_out = BUILD_DATA1(qc5_out.plink, ch_covariates)

    // ── Stage D: ALLFAM stats, PCA, GRM ─────────────────────────────────────

    allfam_out = BUILD_ALLFAM_QC6(qc5_out.plink, ch_fam_2584dogs)
    allfam_prune_out = PRUNE_ALLFAM(allfam_out.plink)
    allfam_pca_out = PCA_ALLFAM(allfam_out.plink, allfam_prune_out.prune_in)
    grm_out = BUILD_GRM(allfam_out.plink)
    depth_hist_out = DEPTH_HISTOGRAMS(ch_fam_qc3modi, ch_depth_file, allfam_out.plink)

    // ── SIZE: control phenotype, run via mlma alongside F1/F2/F3 ────────────
    // Not part of the original pipeline sources — see process_definitions.nf's header comment.

    size_fam_out = BUILD_SIZE_FAM(data1_out.data1, ch_size_pheno_file)
    size_qc6_out = FILTER_SIZE_QC6(size_fam_out.fam, qc5_out.plink)
    size_covar_out = BUILD_SIZE_COVARIATES(size_qc6_out.plink, data1_out.data1)

    // ── Stage E: per-factor phenotype files (F1/F2/F3) ───────────────────────
    // data6 is consumed as a frozen input — see process_definitions.nf's Stage E header for why
    // it isn't recomputed from data1 + F1/F2/F3_CCD3F.txt here.
    //
    // Restricted + reordered to QC5's own dogID set/order before feeding BUILD_FACTOR_FAM/
    // BUILD_ITEM_FAM — data6 is untouched by Stage C's LowPass-12 fix and still has 3328 rows,
    // while QC5 is now 3316; see FILTER_DATA6_TO_QC5's header for why this is required ahead of
    // the positional --fam substitution in FILTER_FACTOR_QC6/FILTER_ITEM_QC6.

    data6_qc5_out = FILTER_DATA6_TO_QC5(ch_data6, qc5_out.plink)

    ch_factor_configs = channel.fromList(getFactorConfigs())
        .map { cfg -> [cfg.factor_id, cfg.item_cols, cfg.na_threshold] }

    factor_fam_out   = BUILD_FACTOR_FAM(ch_factor_configs, data6_qc5_out.data6)
    factor_qc6_out   = FILTER_FACTOR_QC6(factor_fam_out.fam, qc5_out.plink)
    factor_covar_out = BUILD_FACTOR_COVARIATES(factor_qc6_out.plink, data1_out.data1)

    // ── Stage F: per-item phenotype files (14 items) ─────────────────────────

    ch_item_ids = channel.fromList(getItemIds())

    item_fam_out   = BUILD_ITEM_FAM(ch_item_ids, data6_qc5_out.data6)
    item_qc6_out   = FILTER_ITEM_QC6(item_fam_out.fam, qc5_out.plink)
    item_prune_out = PRUNE_ITEM_LD(item_qc6_out.plink)

    // Two independent per-item channels — must be explicitly joined by item_id rather than
    // passed as two separate process inputs, which would pair rows by arrival order, not identity.
    ch_item_for_eigeninput = item_qc6_out.plink.join(item_prune_out.prune_in)

    item_eigeninput_out = BUILD_ITEM_EIGENINPUT(ch_item_for_eigeninput)
    item_grab_out        = BUILD_ITEM_GRAB(item_qc6_out.plink, ch_data6)

    // ── Stage G: cross-phenotype dog/SNP-count check — barrier over all 17 QC6 datasets ──────

    ch_factor_counts = factor_qc6_out.plink.map { factor_id, bed, bim, fam -> [factor_id, fam, bim] }
    ch_item_counts    = item_qc6_out.plink.map { item_id, bed, bim, fam -> ["item${item_id}".toString(), fam, bim] }

    ch_phenotype_counts = ch_factor_counts.mix(ch_item_counts)
        .toSortedList { a, b -> a[0] <=> b[0] }

    phenotype_counts_out = CHECK_PHENOTYPE_COUNTS(
        ch_phenotype_counts.map { rows -> rows.collect { it[0] } },
        ch_phenotype_counts.map { rows -> rows.collect { it[1] } },
        ch_phenotype_counts.map { rows -> rows.collect { it[2] } }
    )

    publish:
    merged_qc3_plink   = flipscan_out.plink
    flipscan_logs      = flipscan_out.flipscan_logs
    exclude_lists      = flipscan_out.exclude_lists
    exclude_round1_log = exclude1_out.log
    exclude_round1_list = exclude1_out.exclude_list
    king_round1_log    = RUN_KING_ROUND1.out.king_log
    king_round1_files  = RUN_KING_ROUND1.out.king_files
    king_round2_log    = RUN_KING_ROUND2.out.king_log
    king_round2_files  = RUN_KING_ROUND2.out.king_files
    merged_qc5_plink   = qc5_out.plink
    data1              = data1_out.data1
    allfam_plink       = allfam_out.plink
    allfam_prune_in    = allfam_prune_out.prune_in
    allfam_prune_out   = allfam_prune_out.prune_out
    allfam_pca_plink   = allfam_pca_out.plink
    allfam_eigenvec    = allfam_pca_out.eigenvec
    allfam_eigenval    = allfam_pca_out.eigenval
    grm_bin            = grm_out.grm_bin
    grm_n_bin          = grm_out.grm_n_bin
    grm_id             = grm_out.grm_id
    depth_histograms_log  = depth_hist_out.log
    depth_histograms_plot = depth_hist_out.plot

    size_fam      = size_fam_out.fam
    size_fam_log  = size_fam_out.log
    size_qc6_plink = size_qc6_out.plink
    size_phen     = size_qc6_out.phen
    size_qcovar   = size_covar_out.qcovar
    size_covar    = size_covar_out.covar
    size_covariates_log = size_covar_out.log

    factor_fam            = factor_fam_out.fam.map { _factor_id, file -> file }
    factor_fam_log        = factor_fam_out.log
    factor_qc6_plink      = factor_qc6_out.plink.map { _factor_id, bed, bim, fam -> [bed, bim, fam] }
    factor_phen           = factor_qc6_out.phen
    factor_qcovar         = factor_covar_out.qcovar
    factor_covar          = factor_covar_out.covar
    factor_covariates_log = factor_covar_out.log
    factor_na_counts_plot = factor_fam_out.na_counts_plot

    item_fam           = item_fam_out.fam.map { _item_id, file -> file }
    item_fam_log       = item_fam_out.log
    item_qc6_plink     = item_qc6_out.plink.map { _item_id, bed, bim, fam -> [bed, bim, fam] }
    item_prune_in      = item_prune_out.prune_in.map { _item_id, file -> file }
    item_prune_out     = item_prune_out.prune_out
    item_eigeninput    = item_eigeninput_out.plink.map { _item_id, bed, bim, fam -> [bed, bim, fam] }
    item_grab          = item_grab_out.grab
    item_grab_log      = item_grab_out.log

    phenotype_counts_log = phenotype_counts_out.log

}

// Output publishing mode is symlink by default. This saves space over copying.
output {
    merged_qc3_plink {
        path 'plink'
    }
    flipscan_logs {
        path 'diagnostics/flipscan'
    }
    exclude_lists {
        path 'diagnostics/flipscan'
    }
    exclude_round1_log {
        path 'diagnostics/exclusion'
    }
    exclude_round1_list {
        path 'diagnostics/exclusion'
    }
    king_round1_log {
        path 'diagnostics/king'
    }
    king_round1_files {
        path 'diagnostics/king'
    }
    king_round2_log {
        path 'diagnostics/king'
    }
    king_round2_files {
        path 'diagnostics/king'
    }
    merged_qc5_plink {
        path 'plink'
    }
    data1 {
        path 'phenotypes'
    }
    allfam_plink {
        path 'plink'
    }
    allfam_prune_in {
        path 'diagnostics/allfam_pca'
    }
    allfam_prune_out {
        path 'diagnostics/allfam_pca'
    }
    allfam_pca_plink {
        path 'plink'
    }
    allfam_eigenvec {
        path 'population_structure'
    }
    allfam_eigenval {
        path 'population_structure'
    }
    grm_bin {
        path 'grm'
    }
    grm_n_bin {
        path 'grm'
    }
    grm_id {
        path 'grm'
    }
    depth_histograms_log {
        path 'diagnostics/depth'
    }
    depth_histograms_plot {
        path 'diagnostics/depth'
    }
    size_fam {
        path 'plink'
    }
    size_fam_log {
        path 'diagnostics/phenotypes'
    }
    size_qc6_plink {
        path 'plink'
    }
    size_phen {
        path 'phenotypes'
    }
    size_qcovar {
        path 'phenotypes'
    }
    size_covar {
        path 'phenotypes'
    }
    size_covariates_log {
        path 'diagnostics/phenotypes'
    }
    factor_fam {
        path 'plink'
    }
    factor_fam_log {
        path 'diagnostics/phenotypes'
    }
    factor_qc6_plink {
        path 'plink'
    }
    factor_phen {
        path 'phenotypes'
    }
    factor_qcovar {
        path 'phenotypes'
    }
    factor_covar {
        path 'phenotypes'
    }
    factor_covariates_log {
        path 'diagnostics/phenotypes'
    }
    factor_na_counts_plot {
        path 'diagnostics/phenotypes'
    }
    item_fam {
        path 'plink'
    }
    item_fam_log {
        path 'diagnostics/phenotypes'
    }
    item_qc6_plink {
        path 'plink'
    }
    item_prune_in {
        path 'diagnostics/phenotypes'
    }
    item_prune_out {
        path 'diagnostics/phenotypes'
    }
    item_eigeninput {
        path 'plink'
    }
    item_grab {
        path 'phenotypes'
    }
    item_grab_log {
        path 'diagnostics/phenotypes'
    }
    phenotype_counts_log {
        path 'diagnostics/phenotypes'
    }
}
