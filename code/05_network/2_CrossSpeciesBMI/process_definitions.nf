// ============================================================================
// process_definitions.nf — 05. Network — 2. CrossSpeciesBMI
//
// Source: 05_network/2_CrossSpeciesBMI/2.1_OCD_dogCD_Network_Colocalization_260521.ipynb
//
// This notebook (and 2.2/2.3) import 3 local helper modules — analysis_functions.py,
// plotting_functions.py, updated_netcoloc_functions.py — via `sys.path.append(cwd); from X import
// *`, committed alongside the notebooks in this directory. Every function used here
// (calculate_network_overlap, calculate_expected_overlap, plot_permutation_histogram,
// get_p_from_permutation_results) is transcribed directly from these files.
// analysis_functions.py's load_pcnet() uses a locally overridden interactome UUID (PCNet2.0).
//
// Reuses 1_NetColoc's container (pandas/numpy/scipy/matplotlib/seaborn/networkx already pinned —
// nothing here needs netcoloc/ddot/ndex-dev specifically).
// ============================================================================

process CONSERVED_NETWORK {
    label 'process_single'
    container 'community.wave.seqera.io/library/netcoloc:8ec66694a6aee19e'

    input:
    path(zcomb_z12, arity: '1')
    path(interactome, arity: '1')

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    conserved_network.py ${zcomb_z12} ${interactome} conserved_network_edgelist.tsv fig1e_threshold_counts.tsv
    """

    stub:
    """
    touch conserved_network_edgelist.tsv fig1e_threshold_counts.tsv
    """

    output:
    path("conserved_network_edgelist.tsv"), emit: edgelist
    path("fig1e_threshold_counts.tsv"), emit: fig1e_counts
}


process PLOT_CONSERVED_NETWORK_VENN {
    // Fig. 1e: NPSh-only / colocalized / NPSd-only 2-circle Venn, from CONSERVED_NETWORK's own
    // fig1e_threshold_counts.tsv -- see bin/plot_conserved_network_venn.R's header for how this
    // was confirmed to match the real published panel exactly. Pure R/eulerr, not part of
    // CONSERVED_NETWORK's own Python environment, so this reuses the network_overlap container
    // built for 05_network/4_Plotting (Fig. 2a-b) instead.
    label 'process_single'
    container 'community.wave.seqera.io/library/network_overlap:282efeb679f77cf4'

    input:
    path(fig1e_counts, arity: '1')

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    plot_conserved_network_venn.R ${fig1e_counts} fig1e_venn.pdf
    """

    stub:
    """
    touch fig1e_venn.pdf
    """

    output:
    path("fig1e_venn.pdf"), emit: plot
}


process NETCOLOC_SIZE_PERMUTATION {
    // Observed network overlap size is fully deterministic; the permuted null distribution is NOT
    // (unseeded random.shuffle in the source) — a fresh permutation is generated each run.
    label 'process_medium'
    container 'community.wave.seqera.io/library/netcoloc:8ec66694a6aee19e'

    input:
    path(zcomb_z12, arity: '1')
    path(dog_seed_genes, arity: '1')
    path(human_seed_genes, arity: '1')

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    netcoloc_size_permutation.py ${zcomb_z12} ${dog_seed_genes} ${human_seed_genes} \\
        netcoloc_size_histogram.pdf observed_network_overlap_size.txt \\
        --num-reps 10000
    """

    stub:
    """
    touch netcoloc_size_histogram.pdf
    touch observed_network_overlap_size.txt
    """

    output:
    path("netcoloc_size_histogram.pdf"), emit: histogram
    path("observed_network_overlap_size.txt"), emit: observed
}


process CONTROL_ANALYSIS_PLOT {
    // controlanalyses_260521.csv is an already-deposited, frozen historical result — the 12-
    // combination control loop that produces it uses an unseeded permutation and is not
    // reproducible, so this process reloads the same file rather than re-deriving it.
    label 'process_single'
    container 'community.wave.seqera.io/library/netcoloc:8ec66694a6aee19e'

    input:
    path(control_results, arity: '1')

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    control_analysis_plot.py ${control_results} control_analyses.pdf
    """

    stub:
    """
    touch control_analyses.pdf
    """

    output:
    path("control_analyses.pdf"), emit: plot
}


// ============================================================================
// 2.2_OCD_dogCD_Systems_Map_260521.ipynb — dogCD/CCD-OCD systems map hierarchy + MGD validation.
// See bin/mgd_significant_communities.py's header for why the notebook's own `final_annotations`
// merge is deliberately not reproduced.
// ============================================================================

process MERGE_HIERARCHY {
    label 'process_single'
    container 'community.wave.seqera.io/library/netcoloc:8ec66694a6aee19e'

    input:
    path(hierarchy, arity: '1')
    path(systems_map, arity: '1')

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    merge_hierarchy.py ${hierarchy} ${systems_map} merged_hierarchy.tsv
    """

    stub:
    """
    touch merged_hierarchy.tsv
    """

    output:
    path("merged_hierarchy.tsv"), emit: merged
}


process MGD_SIGNIFICANT_COMMUNITIES {
    // Reloads the already-deposited, frozen MGD enrichment result rather than recomputing it, which
    // would additionally need the full MGI/MPO/ontology machinery this pipeline doesn't otherwise
    // build. 2 of the 6 deposited files (SCH_CCD, DEP_CCD) use a comma decimal separator instead of
    // a period — coerced explicitly in bin/mgd_significant_communities.py, see its own header.
    tag   "${trait_id}"
    label 'process_single'
    container 'community.wave.seqera.io/library/netcoloc:8ec66694a6aee19e'

    input:
    tuple val(trait_id), path(mgd_enrichment)

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mgd_significant_communities.py ${mgd_enrichment} ${trait_id}_comm_results_sign.tsv
    """

    stub:
    """
    touch ${trait_id}_comm_results_sign.tsv
    """

    output:
    tuple val(trait_id), path("${trait_id}_comm_results_sign.tsv"), emit: comm_results_sign
}


process MGD_NETWORK_ENRICHED_TERMS {
    // Checks the manuscript Results sentence "The vast majority of MP terms enriched in the
    // OCD/dogCD combined network (139 of 142) are associated with C185" against the deposited
    // data -- see bin/mgd_network_enriched_terms.py's own header for why the exact 139/142/10
    // isn't reproducible (that tally was a manual count when Fig. 3c was assembled; neither the
    // frozen MGD file nor its source notebook contains code that computes it) and what this script
    // computes instead (141/133 on the real data, using the same sig_5e6 threshold the deposited
    // file itself already carries).
    //
    // OCD/dogCD-only, unlike MGD_SIGNIFICANT_COMMUNITIES above: this needs a hierarchy file (gene
    // membership per community, to find the root and its immediate children) in addition to the
    // MGD enrichment file, and only 3 of the 7 pairings have one deposited
    // (ocd_ccd_hierarchy/CompulsiveNetwork_hierarchy_data_251022.tsv here; SCH_CCD and OCD_ONLY
    // have their own under different filenames, unused elsewhere in this workflow; DEP_ONLY,
    // SCH_ONLY, DEP_CCD, CCD_ONLY have none) -- generalizing this across all 7 the way the
    // tuple/getMgdPairings() pattern above does isn't possible without deposits that don't exist.
    // Only the OCD/dogCD pairing is reported in the manuscript (Fig. 3c) anyway.
    label 'process_single'
    container 'community.wave.seqera.io/library/netcoloc:8ec66694a6aee19e'

    input:
    path(mgd_enrichment, arity: '1')
    path(hierarchy, arity: '1')

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mgd_network_enriched_terms.py ${mgd_enrichment} ${hierarchy} OCD_CCD_network_enriched_mgd_terms.tsv
    """

    stub:
    """
    touch OCD_CCD_network_enriched_mgd_terms.tsv
    """

    output:
    path("OCD_CCD_network_enriched_mgd_terms.tsv"), emit: network_enriched
}
