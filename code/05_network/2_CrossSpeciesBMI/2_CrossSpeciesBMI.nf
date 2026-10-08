// ============================================================================
// main.nf — 05. Network — 2. CrossSpeciesBMI workflow
//
// dogCD/CCD-OCD conserved network, colocalization size vs. a permuted null (Figure 1g), the
// control-comparisons bar plot, the dogCD/CCD-OCD systems map hierarchy, and per-pairing MGD
// significant-community derivation. See process_definitions.nf's header for the notebooks' local
// helper modules and the `final_annotations` gap.
//
// Source: 05_network/2_CrossSpeciesBMI/{2.1_OCD_dogCD_Network_Colocalization_260521.ipynb,
//          2.2_OCD_dogCD_Systems_Map_260521.ipynb}
// ============================================================================

include { validateParameters } from 'plugin/nf-schema'

include { CONSERVED_NETWORK              } from './process_definitions.nf'
include { PLOT_CONSERVED_NETWORK_VENN    } from './process_definitions.nf'
include { NETCOLOC_SIZE_PERMUTATION      } from './process_definitions.nf'
include { CONTROL_ANALYSIS_PLOT          } from './process_definitions.nf'
include { MERGE_HIERARCHY                } from './process_definitions.nf'
include { MGD_SIGNIFICANT_COMMUNITIES    } from './process_definitions.nf'
include { MGD_NETWORK_ENRICHED_TERMS     } from './process_definitions.nf'

// All 7 systems-map pairings, each with a real, correctly-formatted deposited MGD enrichment
// result — see process_definitions.nf's header for why this is generalized across all of them
// rather than just the OCD-CCD pairing 2.2's own cells literally show. OCD_ONLY's own deposit
// originally had unrelated LDSC annotation content, not MGD enrichment results, in a different
// export format than its siblings; reformat_ocd_only_mgd_results.py normalizes it to match, run
// once to produce the OCD_ONLY_hierarchy_full_MGD_enrichment_results.tsv referenced here.
def getMgdPairings() {
    [
        ['OCD_CCD', '251023_KT_hierarchy_full_MGD_enrichment_results.tsv'],
        ['CCD_ONLY', 'CCD_ONLY_hierarchy_full_MGD_enrichment_results.tsv'],
        ['OCD_ONLY', 'OCD_ONLY_hierarchy_full_MGD_enrichment_results.tsv'],
        ['DEP_ONLY', 'DEP_ONLY_hierarchy_full_MGD_enrichment_results.tsv'],
        ['SCH_ONLY', 'SCH_ONLY_hierarchy_full_MGD_enrichment_results.tsv'],
        ['DEP_CCD', 'DEP_CCD_hierarchy_full_MGD_enrichment_results.tsv'],
        ['SCH_CCD', 'SCH_CCD_hierarchy_full_MGD_enrichment_results.tsv'],
    ]
}

workflow {

    main:
    validateParameters()

    ch_zcomb_z12       = channel.fromPath(params.ocd_ccd_zcomb_z12, checkIfExists: true)
    ch_interactome     = channel.fromPath(params.interactome, checkIfExists: true)
    ch_dog_seed_genes  = channel.fromPath(params.dog_seed_genes, checkIfExists: true)
    ch_human_seed_genes = channel.fromPath(params.human_seed_genes, checkIfExists: true)
    ch_control_results = channel.fromPath(params.control_results, checkIfExists: true)
    ch_systems_map     = channel.fromPath(params.ocd_ccd_systems_map, checkIfExists: true)
    // .collect(): reused below by both MERGE_HIERARCHY and MGD_NETWORK_ENRICHED_TERMS -- a queue
    // channel can only be consumed by one process, so this needs to be a value channel to feed two.
    ch_hierarchy       = channel.fromPath(params.ocd_ccd_hierarchy, checkIfExists: true).collect()

    conserved_out  = CONSERVED_NETWORK(ch_zcomb_z12, ch_interactome)
    venn_out       = PLOT_CONSERVED_NETWORK_VENN(conserved_out.fig1e_counts)
    permutation_out = NETCOLOC_SIZE_PERMUTATION(ch_zcomb_z12, ch_dog_seed_genes, ch_human_seed_genes)
    control_out    = CONTROL_ANALYSIS_PLOT(ch_control_results)

    merged_out = MERGE_HIERARCHY(ch_hierarchy, ch_systems_map)

    ch_mgd_pairings = channel.fromList(getMgdPairings())
        .map { id, fname -> [id, file("${params.mgd_enrichment_dir}/${fname}", checkIfExists: true)] }
    mgd_out = MGD_SIGNIFICANT_COMMUNITIES(ch_mgd_pairings)

    // OCD/dogCD-only (see MGD_NETWORK_ENRICHED_TERMS's own header for why) -- reuses the same
    // filename getMgdPairings() already names for 'OCD_CCD', and the same ch_hierarchy already
    // built above for MERGE_HIERARCHY.
    ch_ocd_ccd_mgd_enrichment = channel.fromPath(
        "${params.mgd_enrichment_dir}/251023_KT_hierarchy_full_MGD_enrichment_results.tsv", checkIfExists: true
    )
    network_enriched_out = MGD_NETWORK_ENRICHED_TERMS(ch_ocd_ccd_mgd_enrichment, ch_hierarchy)

    publish:
    conserved_network_edgelist = conserved_out.edgelist
    fig1e_threshold_counts     = conserved_out.fig1e_counts
    fig1e_venn                 = venn_out.plot
    netcoloc_size_histogram    = permutation_out.histogram
    observed_network_overlap_size = permutation_out.observed
    control_analyses_plot      = control_out.plot
    merged_hierarchy           = merged_out.merged
    comm_results_sign          = mgd_out.comm_results_sign.map { _id, file -> file }
    network_enriched_mgd_terms = network_enriched_out.network_enriched

}

// Output publishing mode is symlink by default. This saves space over copying.
output {
    conserved_network_edgelist {
        path 'network'
    }
    fig1e_threshold_counts {
        path 'network'
    }
    fig1e_venn {
        path 'network'
    }
    netcoloc_size_histogram {
        path 'network'
    }
    observed_network_overlap_size {
        path 'network'
    }
    control_analyses_plot {
        path 'controls'
    }
    merged_hierarchy {
        path 'systems_map'
    }
    comm_results_sign {
        path 'mgd_enrichment'
    }
    network_enriched_mgd_terms {
        path 'mgd_enrichment'
    }
}
