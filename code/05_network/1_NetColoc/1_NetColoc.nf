// ============================================================================
// main.nf — 05. Network — 1. NetColoc workflow
//
// Fetches PCNet2.0 once, builds the (fully deterministic) propagation matrices once, then
// computes network-propagation z-scores for each of the 4 distinct dogCD/human-trait gene sets —
// see process_definitions.nf's header for why this differs in shape from the original notebooks'
// own redundant per-cross-species-pairing recomputation.
//
// Source: 05_network/1_NetColoc/{1.1_OCD_dogCD_NetColoc_analysis_260521.ipynb,
//          1.2_single-species_OCD_NetColoc_analysis_260521.ipynb}
// ============================================================================

include { validateParameters } from 'plugin/nf-schema'

include { FETCH_PCNET2          } from './process_definitions.nf'
include { BUILD_HEATS_MATRIX    } from './process_definitions.nf'
include { COMPUTE_ZSCORES       } from './process_definitions.nf'
include { COMBINE_ZSCORES       } from './process_definitions.nf'
include { PLOT_NETCOLOC_SCATTER } from './process_definitions.nf'

// The 4 distinct gene sets actually z-scored across both notebooks (OCD/DEP/SCH from 1.2's
// single-species table, dogCD/CCD common to all of 1.1's cross-species pairings) — filenames
// match exactly what's read in notebook cells 6/7 (1.1) and the equivalent cells in 1.2.
def getTraits() {
    [
        ['OCD', 'OCDgenes_v4.txt'],
        ['CCD', 'CCDgenes_Oct25.txt'],
        ['DEP', 'human_DEP_noMHC_seed_genes_header.txt'],
        ['SCH', 'human_schiz_seed_genes_header.txt'],
    ]
}

workflow {

    main:
    validateParameters()

    ch_uuid = channel.value(params.pcnet_uuid)
    fetch_out = FETCH_PCNET2(ch_uuid)

    matrix_out = BUILD_HEATS_MATRIX(fetch_out.interactome)

    ch_seed_genes = channel.fromList(getTraits())
        .map { id, fname -> [id, file("${params.seed_genes_dir}/${fname}", checkIfExists: true)] }

    // .collect() converts these single-emission process-output channels into broadcastable value
    // channels — see 04_gwas/3_finemap_susie's own established rationale for preferring this over
    // .first() (which would silently drop the rest if either ever unexpectedly emitted more than
    // once, masking a real bug, rather than surfacing it).
    ch_interactome_bc = fetch_out.interactome.collect()
    ch_indiv_heats_bc = matrix_out.indiv_heats_matrix.collect()

    zscores_out = COMPUTE_ZSCORES(ch_seed_genes, ch_interactome_bc, ch_indiv_heats_bc)

    // dogCD/CCD is the constant "r" side of every deposited *_zcomb_z12* file (cell 24's
    // full_join(ccd, ocd, ...) plus the two undocumented-but-deposited DEP/SCH pairings) — combine
    // it against each of the other 3 human traits (OCD, DEP, SCH).
    ch_ccd_zscore_bc = zscores_out.z_scores
        .filter { id, _f -> id == 'CCD' }
        .map { _id, f -> f }
        .collect()

    ch_human_trait_zscores = zscores_out.z_scores
        .filter { id, _f -> id != 'CCD' }

    zcomb_out = COMBINE_ZSCORES(ch_human_trait_zscores, ch_ccd_zscore_bc)

    // Fig. 1f only needs the OCD/dogCD pairing's own combined z-scores (not DEP/SCH's) — same
    // single pairing CONSERVED_NETWORK/etc. in 2_CrossSpeciesBMI are scoped to.
    ch_ocd_ccd_zcomb = zcomb_out.zcomb
        .filter { id, _f -> id == 'OCD' }
        .map { _id, f -> f }

    ch_toga_file = channel.fromPath(params.toga_file, checkIfExists: true)

    ch_ocd_seed_genes = ch_seed_genes
        .filter { id, _f -> id == 'OCD' }
        .map { _id, f -> f }
    ch_ccd_seed_genes = ch_seed_genes
        .filter { id, _f -> id == 'CCD' }
        .map { _id, f -> f }

    scatter_out = PLOT_NETCOLOC_SCATTER(ch_ocd_ccd_zcomb, ch_toga_file, ch_ocd_seed_genes, ch_ccd_seed_genes)

    publish:
    interactome        = fetch_out.interactome
    w_prime             = matrix_out.w_prime
    indiv_heats_matrix  = matrix_out.indiv_heats_matrix
    z_scores            = zscores_out.z_scores.map { _id, file -> file }
    zcomb_z12           = zcomb_out.zcomb.map { _id, file -> file }
    fig1f_combined_table = scatter_out.combined_table
    fig1f_boxplot        = scatter_out.boxplot
    fig1f_source_data    = scatter_out.source_data
    fig1f_scatter        = scatter_out.scatter

}

// Output publishing mode is symlink by default. This saves space over copying.
output {
    interactome {
        path 'interactome'
    }
    w_prime {
        path 'matrices'
    }
    indiv_heats_matrix {
        path 'matrices'
    }
    z_scores {
        path 'z_scores'
    }
    zcomb_z12 {
        path 'zcomb_z12'
    }
    fig1f_combined_table {
        path 'fig1f'
    }
    fig1f_boxplot {
        path 'fig1f'
    }
    fig1f_source_data {
        path 'fig1f'
    }
    fig1f_scatter {
        path 'fig1f'
    }
}
