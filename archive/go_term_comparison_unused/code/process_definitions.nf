// ============================================================================
// process_definitions.nf — 05. Network — 4. Plotting
//
// Source: 05_network/4_Plotting/GO-term_comparison.ipynb
//
// Reuses 1_NetColoc's container (pandas/numpy/matplotlib are already pinned there) rather than
// building a new one for a 3-cell plotting notebook with no netcoloc-specific dependency.
// ============================================================================

process GO_TERM_COMPARISON {
    // One task per trait (OCD, DEP, SCH) — see bin/go_term_comparison.py's header for the
    // xlim/ylim divergence between OCD (dynamic) and DEP/SCH (hardcoded) reproduced literally.
    // Output filenames are passed in explicitly (not derived from trait_id) because the original
    // notebook's own naming is inconsistent across traits (e.g. "GO_enrichment_CCDOCD_OCD_..." vs
    // "GO_enrichment_CCDDEP-DEP_..." — underscore vs hyphen) — reproduced literally rather than
    // normalized.
    tag   "${trait_id}"
    label 'process_single'
    container 'community.wave.seqera.io/library/netcoloc:8ec66694a6aee19e'

    input:
    tuple val(trait_id), path(run1_go_enrichment), path(run2_go_enrichment), val(xlim), val(ylim), val(run1_only_name), val(full_table_name), val(plot_name)

    when:
    task.ext.when == null || task.ext.when

    script:
    def xlim_opt = xlim ? "--xlim ${xlim[0]} ${xlim[1]}" : ''
    def ylim_opt = ylim ? "--ylim ${ylim[0]} ${ylim[1]}" : ''
    """
    go_term_comparison.py ${run1_go_enrichment} ${run2_go_enrichment} ${trait_id} \\
        ${run1_only_name} \\
        ${full_table_name} \\
        ${plot_name} \\
        ${xlim_opt} ${ylim_opt}
    """

    stub:
    """
    touch ${run1_only_name}
    touch ${full_table_name}
    touch ${plot_name}
    """

    output:
    tuple val(trait_id), path("${run1_only_name}"), emit: run1_only_terms
    tuple val(trait_id), path("${full_table_name}"), emit: full_table
    tuple val(trait_id), path("${plot_name}"), emit: plot
}
