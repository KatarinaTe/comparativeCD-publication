// ============================================================================
// main.nf — 1. Survey data workflow (EFA/IRT)
//
// 14-item dogCD survey responses -> EFA diagnostics -> per-factor IRT fits (F1/F2/F3) ->
// per-dog factor scores, consumed by 3_Merging_filtering.
//
// Source: 02_data_processing/1_Survey_data/Darwin_dog-CD_3-factorsolution.R
// ============================================================================

include { validateParameters } from 'plugin/nf-schema'

include { RUN_EFA                    } from './process_definitions.nf'
include { RUN_IRT_FACTOR             } from './process_definitions.nf'
include { PLOT_EXTRA_ITEM_SE         } from './process_definitions.nf'
include { PLOT_FACTOR_DISTRIBUTIONS  } from './process_definitions.nf'

// Fixed per-factor configuration — a pre-decided grouping/parameterization from the original
// analysis, not a runtime input. See bin/run_irt_factor.R's header for what each field preserves.
def getFactorConfigs() {
    [
        [
            factor_id:            'F1',
            item_cols:            '7,153,154,155',
            item_weights:         '0.81,4.86,25.77,1.44',
            m2_type:              'M2*',
            score_method:         'MAP',
            extra_item_plot_col:  '7',
            trace_theta_range:    '-3,8',
            testinfo_theta_range: '-3,3',
            iteminfo_theta_range: '-3,3',
            iteminfo_palette:     'Set3',
        ],
        [
            factor_id:            'F2',
            item_cols:            '93,95,150',
            item_weights:         '1.01,2.69,1.58',
            m2_type:              'C2',
            score_method:         'EAP',
            extra_item_plot_col:  '',
            trace_theta_range:    '-3,8',
            testinfo_theta_range: '-4,4',
            iteminfo_theta_range: '-4,6',
            iteminfo_palette:     'Set2',
        ],
        [
            factor_id:            'F3',
            item_cols:            '145,146,147,148,149,151,152',
            item_weights:         '0.83,1.12,1.21,1.07,1.12,1.10,1.12',
            m2_type:              'M2*',
            score_method:         'MAP',
            extra_item_plot_col:  '',
            trace_theta_range:    '-3,12',
            testinfo_theta_range: '-5,13',
            iteminfo_theta_range: '-5,13',
            iteminfo_palette:     'Set3',
        ],
    ]
}

workflow {

    main:
    validateParameters()

    ch_response_df    = channel.value(file(params.response_df,           checkIfExists: true))
    ch_questions_csv  = channel.value(file(params.survey_questions_csv,  checkIfExists: true))

    // ── EFA diagnostics and fit — reproduced in full, feeds nothing downstream ──────────

    efa_out = RUN_EFA(ch_response_df, ch_questions_csv)

    // ── Per-factor IRT fits ──────────────────────────────────────────────────────────────

    factor_configs = getFactorConfigs()

    ch_factor_configs = channel.fromList(factor_configs)
        .map { cfg -> [cfg.factor_id, cfg.item_cols, cfg.item_weights, cfg.m2_type, cfg.score_method,
                       cfg.trace_theta_range, cfg.testinfo_theta_range, cfg.iteminfo_theta_range, cfg.iteminfo_palette] }

    irt_out = RUN_IRT_FACTOR(ch_factor_configs, ch_response_df)

    // ── F1's one extra diagnostic plot — invoked only for the factor(s) configured for it ──

    ch_extra_item_plot_cols = channel.fromList(factor_configs)
        .map { cfg -> [cfg.factor_id, cfg.extra_item_plot_col] }
        .filter { _factor_id, extra_item_plot_col -> extra_item_plot_col != '' }

    ch_for_extra_plot = irt_out.scores.join(ch_extra_item_plot_cols)

    extra_plot_out = PLOT_EXTRA_ITEM_SE(ch_for_extra_plot)

    // Order F1/F2/F3 deterministically before the final barrier — factor completion order on a
    // cluster is not guaranteed to be F1,F2,F3.
    ch_ordered_scores = irt_out.scores
        .toSortedList { a, b -> a[0] <=> b[0] }
        .map { sorted -> sorted.collect { it[1] } }

    // ── Final distribution plots — needs all three factors together ────────────────────

    dist_out = PLOT_FACTOR_DISTRIBUTIONS(ch_ordered_scores)

    ch_all_versions = channel.topic("versions")
        .unique()
        .map { process, tool, version ->
            [ process.tokenize(":").last(), "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }
        .collectFile(name: 'versions.yml', sort: true, newLine: true)

    publish:
    factor_scores            = irt_out.scores.map { _factor_id, file -> file }
    efa_log                  = efa_out.log
    efa_parallel_plot        = efa_out.parallel_analysis_plot
    efa_nscree_plot          = efa_out.nscree_plot
    efa_uniqueness_plot      = efa_out.uniqueness_plot
    irt_log                  = irt_out.log
    irt_trace_plot           = irt_out.trace_plot
    irt_iteminfo_facet_plot  = irt_out.iteminfo_facet_plot
    irt_testinfo_plot        = irt_out.testinfo_plot
    irt_iteminfo_plot        = irt_out.iteminfo_plot
    irt_eap_map_plot         = irt_out.eap_map_plot
    irt_na_count_vs_se_plot     = irt_out.na_count_vs_se_plot
    irt_na_count_vs_score_plot  = irt_out.na_count_vs_score_plot
    irt_score_vs_se_plot        = irt_out.score_vs_se_plot
    extra_item_plot          = extra_plot_out.plot
    distribution_log         = dist_out.log
    distribution_f1_plot     = dist_out.f1_plot
    distribution_f2_plot     = dist_out.f2_plot
    distribution_f3_plot     = dist_out.f3_plot
    versions                 = ch_all_versions

}

// Output publishing mode is symlink by default. This saves space over copying.
output {
    factor_scores {
        path 'scores'
    }
    efa_log {
        path 'efa/diagnostics'
    }
    efa_parallel_plot {
        path 'efa/plots'
    }
    efa_nscree_plot {
        path 'efa/plots'
    }
    efa_uniqueness_plot {
        path 'efa/plots'
    }
    irt_log {
        path 'irt/diagnostics'
    }
    irt_trace_plot {
        path 'irt/plots'
    }
    irt_iteminfo_facet_plot {
        path 'irt/plots'
    }
    irt_testinfo_plot {
        path 'irt/plots'
    }
    irt_iteminfo_plot {
        path 'irt/plots'
    }
    irt_eap_map_plot {
        path 'irt/plots'
    }
    irt_na_count_vs_se_plot {
        path 'irt/plots'
    }
    irt_na_count_vs_score_plot {
        path 'irt/plots'
    }
    irt_score_vs_se_plot {
        path 'irt/plots'
    }
    extra_item_plot {
        path 'irt/plots'
    }
    distribution_log {
        path 'distributions/diagnostics'
    }
    distribution_f1_plot {
        path 'distributions/plots'
    }
    distribution_f2_plot {
        path 'distributions/plots'
    }
    distribution_f3_plot {
        path 'distributions/plots'
    }
    versions {
        path 'pipeline_info'
    }
}
