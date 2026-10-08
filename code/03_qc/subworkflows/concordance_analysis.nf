include { ISEC_COMMON_SITES     } from '../process_definitions.nf'
include { RESTRICT_TO_SITES_VCF } from '../process_definitions.nf'
include { UPDATE_VCF_DICT       } from '../process_definitions.nf'
include { GENOTYPE_CONCORDANCE  } from '../process_definitions.nf'

// Comparison tail shared by all 12 concordance comparisons in this pipeline (study 1, study 2,
// and one per downsampling fraction), called once over a channel mixing all of them together,
// tagged by label. The call-side and truth-side VCF are likewise mixed into one channel (tagged
// "<label>_call" / "<label>_truth") before RESTRICT_TO_SITES_VCF/UPDATE_VCF_DICT.
workflow CONCORDANCE_ANALYSIS {
    take:
    comparisons // channel: [label, call_vcf, call_idx, truth_vcf, truth_idx] — one row per comparison
    pairs       // channel: [label, pair_id, call_sample, truth_sample] — many rows per label
    ref_dict    // channel: path — single emission

    main:
    sites_out = ISEC_COMMON_SITES(comparisons)

    ch_call_with_sites = comparisons
        .map { label, call_vcf, call_idx, _truth_vcf, _truth_idx -> [label, call_vcf, call_idx] }
        .join(sites_out.sites)
        .map { label, vcf, idx, sites -> ["${label}_call", vcf, idx, sites] }
    ch_truth_with_sites = comparisons
        .map { label, _call_vcf, _call_idx, truth_vcf, truth_idx -> [label, truth_vcf, truth_idx] }
        .join(sites_out.sites)
        .map { label, vcf, idx, sites -> ["${label}_truth", vcf, idx, sites] }

    restricted = RESTRICT_TO_SITES_VCF(ch_call_with_sites.mix(ch_truth_with_sites))
    dict_out   = UPDATE_VCF_DICT(restricted.vcf, ref_dict)

    // Strip the _call/_truth suffix added above so results can be joined back to `pairs` by
    // the original comparison label. Anchored to end-of-string — safe as long as no label
    // itself ends in "_call"/"_truth" (none of this pipeline's labels do).
    ch_call_final  = dict_out.vcf
        .filter { label, _vcf, _idx -> label.endsWith('_call') }
        .map    { label, vcf, idx -> [label.replaceAll(/_call$/, ''), vcf, idx] }
    ch_truth_final = dict_out.vcf
        .filter { label, _vcf, _idx -> label.endsWith('_truth') }
        .map    { label, vcf, idx -> [label.replaceAll(/_truth$/, ''), vcf, idx] }

    concordance_in = pairs
        .combine(ch_call_final,  by: 0)
        .combine(ch_truth_final, by: 0)
        .map { _label, pair_id, call_sample, truth_sample, call_vcf, call_idx, truth_vcf, truth_idx ->
            [pair_id, call_sample, truth_sample, call_vcf, call_idx, truth_vcf, truth_idx]
        }

    concordance_out = GENOTYPE_CONCORDANCE(concordance_in)

    emit:
    summary_metrics     = concordance_out.summary_metrics
    contingency_metrics = concordance_out.contingency_metrics
    detail_metrics      = concordance_out.detail_metrics
}
