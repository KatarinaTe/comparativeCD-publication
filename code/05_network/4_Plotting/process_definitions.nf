// ============================================================================
// process_definitions.nf — 05. Network — 4. Plotting (Fig. 2a-b)
//
// Source: 05_network/4_Plotting/Fig2_network-overlap.R -- see
// archive/fig2b_euler_diagram_superseded/ for the original Panel A/B conversion story, including
// the real bug found and fixed 2026-09-22 (Panel B was reading the wrong -- cross-species instead
// of human-only -- Depression/Schizophrenia gene lists).
//
// Panel B replaced 2026-10-01: the stacked-bar + Euler diagram is superseded by area-proportional
// 3-circle Venn diagrams (bin/plot_venn_hierarchy.R) -- see that script's header and this
// directory's README for the conversion/verification story.
// ============================================================================

process PLOT_VENN_HIERARCHY {
    // Fig. 2b -- area-proportional Venn diagrams of in_hierarchy genes (dogCD / OCD / OCD-dogCD
    // vs. Depression/Schizophrenia). Reuses go_semantic_plots (tidyverse + cowplot already
    // pinned there) -- no eulerr dependency, so no new container needed.
    label 'process_single'
    container 'community.wave.seqera.io/library/go_semantic_plots:69629d5c597a400b'

    input:
    path(gene_level_network_basis_table, arity: '1')

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    plot_venn_hierarchy.R ${gene_level_network_basis_table}
    """

    stub:
    """
    touch venn_hierarchy_dogCD.pdf venn_hierarchy_dogCD.png
    touch venn_hierarchy_OCD.pdf venn_hierarchy_OCD.png
    touch venn_hierarchy_OCD_dogCD.pdf venn_hierarchy_OCD_dogCD.png
    """

    output:
    path("venn_hierarchy_dogCD.pdf"),     emit: dogcd_pdf
    path("venn_hierarchy_dogCD.png"),     emit: dogcd_png
    path("venn_hierarchy_OCD.pdf"),       emit: ocd_pdf
    path("venn_hierarchy_OCD.png"),       emit: ocd_png
    path("venn_hierarchy_OCD_dogCD.pdf"), emit: ocd_dogcd_pdf
    path("venn_hierarchy_OCD_dogCD.png"), emit: ocd_dogcd_png
}


process ADD_VENN_REGION_LABELS {
    // Adds a human-readable "which Fig. 2b Venn region" column per comparison (dogCD/OCD/
    // OCD-dogCD vs. Depression/Schizophrenia) to the gene-level table (Supplementary Table 8's
    // source data) -- e.g. "Depression + Schizophrenia", "OCD only". Reuses the exact same region
    // definitions PLOT_VENN_HIERARCHY draws from (same in_hierarchy_* columns, same 3
    // comparisons), so the new columns are guaranteed consistent with the published Venns by
    // construction, not a separate reimplementation that could drift. Verified to reproduce the
    // same 7 region counts per comparison as the published figure -- see
    // bin/add_venn_region_labels.R's own header.
    label 'process_single'
    container 'community.wave.seqera.io/library/go_semantic_plots:69629d5c597a400b'

    input:
    path(gene_level_network_basis_table, arity: '1')

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    add_venn_region_labels.R ${gene_level_network_basis_table} gene_level_network_basis_table_with_venn_regions.csv
    """

    stub:
    """
    touch gene_level_network_basis_table_with_venn_regions.csv
    """

    output:
    path("gene_level_network_basis_table_with_venn_regions.csv"), emit: table
}


process BUILD_NETWORK_GENE_COUNTS {
    // Fig. 2a's underlying data -- live-computed replacement for the hardcoded (captured, total)
    // literals plot_network_overlap.R's Panel A used to draw. See bin/build_network_gene_counts.py's
    // header for the real bug this closes (Depression/dogCD's total was coded assuming the wrong
    // z_comb threshold; a second bug, a hand-tracked "8 dog seeds" copied across all three
    // cross-species rows, is also fixed by computing this live instead).
    //
    // Reuses 1_NetColoc's netcoloc container (pandas already pinned there) -- pure data
    // wrangling, no netcoloc-specific functions actually called.
    label 'process_single'
    container 'community.wave.seqera.io/library/netcoloc:8ec66694a6aee19e'

    input:
    path(dogcd_z, arity: '1')
    path(dogcd_seeds, arity: '1')
    path(dogcd_hier, arity: '1')
    path(ocd_z, arity: '1')
    path(ocd_seeds, arity: '1')
    path(ocd_hier, arity: '1')
    path(dep_z, arity: '1')
    path(dep_seeds, arity: '1')
    path(dep_hier, arity: '1')
    path(sch_z, arity: '1')
    path(sch_seeds, arity: '1')
    path(sch_hier, arity: '1')
    path(ocd_ccd_zcomb, arity: '1')
    path(ocd_ccd_hier, arity: '1')
    path(dep_ccd_zcomb, arity: '1')
    path(dep_ccd_hier, arity: '1')
    path(sch_ccd_zcomb, arity: '1')
    path(sch_ccd_hier, arity: '1')

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    build_network_gene_counts.py \\
        --dogcd-z ${dogcd_z} --dogcd-seeds ${dogcd_seeds} --dogcd-seeds-header --dogcd-hier ${dogcd_hier} \\
        --ocd-z ${ocd_z} --ocd-seeds ${ocd_seeds} --ocd-seeds-header --ocd-hier ${ocd_hier} \\
        --dep-z ${dep_z} --dep-seeds ${dep_seeds} --dep-hier ${dep_hier} \\
        --sch-z ${sch_z} --sch-seeds ${sch_seeds} --sch-hier ${sch_hier} \\
        --ocd-ccd-zcomb ${ocd_ccd_zcomb} --ocd-ccd-hier ${ocd_ccd_hier} \\
        --dep-ccd-zcomb ${dep_ccd_zcomb} --dep-ccd-hier ${dep_ccd_hier} \\
        --sch-ccd-zcomb ${sch_ccd_zcomb} --sch-ccd-hier ${sch_ccd_hier} \\
        network_gene_counts.csv
    """

    stub:
    """
    touch network_gene_counts.csv
    """

    output:
    path("network_gene_counts.csv"), emit: counts
}


process PLOT_FIG2A {
    // Reuses go_semantic_plots (tidyverse + cowplot already pinned there, same container
    // 07e-power-analysis reuses) -- no new build.
    label 'process_single'
    container 'community.wave.seqera.io/library/go_semantic_plots:69629d5c597a400b'

    input:
    path(gene_counts, arity: '1')

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    plot_fig2a.R ${gene_counts} fig2a.pdf
    """

    stub:
    """
    touch fig2a.pdf
    """

    output:
    path("fig2a.pdf"), emit: plot
}
