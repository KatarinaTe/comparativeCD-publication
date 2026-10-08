// ============================================================================
// main.nf — 05. Network — 4. Plotting workflow
//
// GO-enrichment "run1 (single-species-only) vs run2 (dogCD/CCD cross-species)" difference volcano
// plots, for each of the 3 human traits (OCD, DEP, SCH).
//
// Source: 05_network/4_Plotting/GO-term_comparison.ipynb
// ============================================================================

include { validateParameters } from 'plugin/nf-schema'

include { GO_TERM_COMPARISON } from './process_definitions.nf'

// One entry per trait: [trait_id, run1 (ONLY) file, run2 (CCD-paired) file, xlim, ylim,
// run1-only-terms filename, full-table filename, plot filename] — xlim/ylim are null for OCD
// (computed dynamically, matching cell 1) and fixed for DEP/SCH (matching cells 2/3's hardcoded
// values). Filenames match the notebook's own literal (and inconsistent) naming exactly.
def getComparisons() {
    [
        ['OCD', 'OCD_ONLY_hierachy_full_GO_enrichment.tsv', 'Compulsive_hierachy_full_GO_enrichment_251023.tsv', null, null,
         'OCD_run1_only_terms.txt', 'CCD-OCD_OCD_GO_difference_volcano_all_terms.txt', 'GO_enrichment_CCDOCD_OCD_volcano_color.pdf'],
        ['DEP', 'DEP_ONLY_hierarchy_full_GO_enrichment.tsv', 'DEP_CCD_hierachy_full_GO_enrichment.tsv', [-20, 70], [0, 140],
         'DEP_run1_only_terms.txt', 'CCD-DEP_DEP_GO_difference_volcano_all_terms.txt', 'GO_enrichment_CCDDEP-DEP_volcano_color.pdf'],
        ['SCH', 'SCH_ONLY_hierachy_full_GO_enrichment.tsv', 'SCH_CCD_hierachy_full_GO_enrichment.tsv', [-20, 70], [0, 140],
         'SCH_run1_only_terms.txt', 'CCD-SCH_SCH_GO_difference_volcano_all_terms.txt', 'GO_enrichment_CCDSCH-SCH_volcano_color.pdf'],
    ]
}

workflow {

    main:
    validateParameters()

    ch_comparisons = channel.fromList(getComparisons())
        .map { trait_id, run1_name, run2_name, xlim, ylim, run1_only_name, full_table_name, plot_name ->
            [
                trait_id,
                file("${params.datadir}/${run1_name}", checkIfExists: true),
                file("${params.datadir}/${run2_name}", checkIfExists: true),
                xlim, ylim, run1_only_name, full_table_name, plot_name,
            ]
        }

    go_out = GO_TERM_COMPARISON(ch_comparisons)

    publish:
    run1_only_terms = go_out.run1_only_terms.map { _id, file -> file }
    full_table       = go_out.full_table.map { _id, file -> file }
    plot             = go_out.plot.map { _id, file -> file }

}

// Output publishing mode is symlink by default. This saves space over copying.
output {
    run1_only_terms {
        path 'plotting'
    }
    full_table {
        path 'plotting'
    }
    plot {
        path 'plotting'
    }
}
