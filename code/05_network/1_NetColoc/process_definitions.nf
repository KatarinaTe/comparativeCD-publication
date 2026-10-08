// ============================================================================
// process_definitions.nf — 05. Network — 1. NetColoc
//
// Source: 05_network/1_NetColoc/{1.1_OCD_dogCD_NetColoc_analysis_260521.ipynb,
//          1.2_single-species_OCD_NetColoc_analysis_260521.ipynb}
//
// The notebooks recompute the same 4 gene sets' z-scores redundantly across separate kernel
// sessions (dogCD's z-scores alone get computed 3x across the 3 cross-species pairings in 1.1,
// plus once more in 1.2) — an artifact of the notebooks' own interactive, cell-by-cell structure.
// This pipeline computes each of the 4 distinct gene sets (OCD, dogCD/CCD, DEP, SCH) exactly
// once; downstream cross-species combination steps (2_CrossSpeciesBMI) reuse these same z-score
// outputs by reference.
//
// Two of these 4 gene-set computations are fully reproducible:
//   - w_prime/individual_heats_matrix (get_normalized_adjacency_matrix/get_individual_heats_matrix):
//     pure deterministic linear algebra, no randomness anywhere. Recomputed live here rather than
//     requiring the deposited 2.8GB precomputed matrices the original README offers as a shortcut
//     — genuinely cheap (a dense O(N^3) inversion, a few minutes on one node with decent BLAS).
//   - calculate_heat_zscores: has a `random_seed=1` DEFAULT parameter, called via
//     `np.random.seed(random_seed)` as the function's very first line. Neither notebook ever
//     overrides it, so every historical z-score run used this same fixed seed — reproduced
//     explicitly here.
//
// New container required: both NetColoc and CrossSpeciesBMI (github.com/ucsd-ccbb/NetColoc,
// github.com/sarah-n-wright/CrossSpeciesBMI) are external Python tools. Built from
// CrossSpeciesBMI's own real, fully-pinned environment.yml — see env/05_network/README.md.
// ============================================================================

process FETCH_PCNET2 {
    // Live NDEx fetch (public, no credentials) — the raw CX export for PCNet2.0 is >1GB, larger
    // than the "19,267 nodes" summary figure suggests, and the connection is intermittently flaky
    // for a transfer this size — retried with backoff inside bin/fetch_pcnet2.py itself, rather
    // than via Nextflow's own task-level retry (which would re-run the whole task just to retry
    // one HTTP call).
    label 'process_medium'
    container 'community.wave.seqera.io/library/netcoloc:8ec66694a6aee19e'

    input:
    val(uuid)

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    fetch_pcnet2.py --uuid ${uuid} interactome.pkl
    """

    stub:
    """
    touch interactome.pkl
    """

    output:
    path("interactome.pkl"), emit: interactome
}


process BUILD_HEATS_MATRIX {
    // The O(N^3) dense matrix inversion — see this file's header for why this is recomputed live
    // rather than requiring the deposited matrices.
    label 'process_high'
    container 'community.wave.seqera.io/library/netcoloc:8ec66694a6aee19e'

    input:
    path(interactome, arity: '1')

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    build_heats_matrix.py ${interactome} w_prime.npy indiv_heats_matrix.npy
    """

    stub:
    """
    touch w_prime.npy
    touch indiv_heats_matrix.npy
    """

    output:
    path("w_prime.npy"), emit: w_prime
    path("indiv_heats_matrix.npy"), emit: indiv_heats_matrix
}


process COMPUTE_ZSCORES {
    // One task per distinct gene set (4 total: OCD, dogCD/CCD, DEP, SCH) — see this file's header
    // for why this doesn't match the notebooks' own redundant per-pairing recomputation.
    tag   "${trait_id}"
    label 'process_high'
    container 'community.wave.seqera.io/library/netcoloc:8ec66694a6aee19e'

    input:
    tuple val(trait_id), path(seed_gene_file)
    path(interactome, arity: '1')
    path(indiv_heats_matrix, arity: '1')

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    compute_zscores.py ${interactome} ${indiv_heats_matrix} ${seed_gene_file} ${trait_id}_z_scores.csv
    """

    stub:
    """
    touch ${trait_id}_z_scores.csv
    """

    output:
    tuple val(trait_id), path("${trait_id}_z_scores.csv"), emit: z_scores
}


process COMBINE_ZSCORES {
    // One task per human trait paired against dogCD/CCD (OCD, DEP, SCH) — see bin/combine_zscores.py's
    // header for the "r"/"h"/"hr" column-naming story.
    tag   "${trait_id}"
    label 'process_single'
    container 'community.wave.seqera.io/library/netcoloc:8ec66694a6aee19e'

    input:
    tuple val(trait_id), path(h_zscores)
    path(r_zscores, arity: '1')

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    combine_zscores.py ${r_zscores} ${h_zscores} ${trait_id}_CCD_zcomb_z12.tsv
    """

    stub:
    """
    touch ${trait_id}_CCD_zcomb_z12.tsv
    """

    output:
    tuple val(trait_id), path("${trait_id}_CCD_zcomb_z12.tsv"), emit: zcomb
}


process PLOT_NETCOLOC_SCATTER {
    // Fig. 1f: dogCD-vs-OCD z-score scatter plot, plus a diagnostic orthology-class boxplot —
    // see bin/plot_netcoloc_zscores.R's header for the full provenance/conversion story. Uses the
    // shared R container (env/shared/gwas_supplement_plots.yml), not netcoloc -- this step is
    // pure R/ggplot2, no NetColoc/Python tooling involved.
    label 'process_single'
    container 'community.wave.seqera.io/library/gwas_supplement_plots:e626727dcb23673b'

    input:
    path(ocd_ccd_zcomb, arity: '1')
    path(toga_file, arity: '1')
    path(human_seed_genes, arity: '1')
    path(dog_seed_genes, arity: '1')

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    plot_netcoloc_zscores.R ${ocd_ccd_zcomb} ${toga_file} ${human_seed_genes} ${dog_seed_genes} \\
        ocd_ccd_tab.tsv orthology_class_boxplot.pdf Figure_1f_source_data.tsv fig1f_scatter.pdf
    """

    stub:
    """
    touch ocd_ccd_tab.tsv orthology_class_boxplot.pdf Figure_1f_source_data.tsv fig1f_scatter.pdf
    """

    output:
    path("ocd_ccd_tab.tsv"), emit: combined_table
    path("orthology_class_boxplot.pdf"), emit: boxplot
    path("Figure_1f_source_data.tsv"), emit: source_data
    path("fig1f_scatter.pdf"), emit: scatter
}
