# data/03_qc

Sample/pair lists for [`code/03_qc`](../../code/03_qc), plus frozen historical output copies from
the original (pre-conversion) concordance runs. The original pre-conversion `README_03.qc.md` was
already superseded and archived before this pass (empty stub found in this location, full content
already in [`archive/legacy_readmes/03_qc_README_03.qc.md`](../../archive/legacy_readmes/03_qc_README_03.qc.md)).

## Live pipeline inputs

| File | Consumed by |
|---|---|
| `samples_to_extract.txt` | `samples_highpass_11` param — Study 1 (HighPass vs. LowPass) |
| `samples_to_extract_down.txt` | `samples_downsampling_10` param — Study 3 (downsampling) |
| `downsample_p_at_1x.csv` | `downsample_p_at_1x` param — Study 3 |
| `concordance_highpass_lowpass_pairs.csv` | `concordance_highpass_lowpass_pairs` param — Study 1 (15 pairs) |
| `concordance_axiom_lowpass_pairs.csv` | `concordance_axiom_lowpass_pairs` param — Study 2 (6 pairs, incl. the deliberate `dogNO_PAIR` negative control) |

## Frozen historical outputs (from the original pre-conversion runs, not live inputs)

`*.v2_concordance.vcf.genotype_concordance_summary_metrics` (Study 1),
`wholedog_gc_*_concordance.vcf.genotype_concordance_contingency_metrics` (Study 2),
`combined_concordance.txt`/`combined_knownsites.depth.txt` (Study 3) — kept for reference/
provenance; the live pipeline regenerates its own equivalents under
`data/results/03_qc/concordance/`.

## Other

- `key_LOW_vs_HIGH.csv`, `DA_downsampling_bam_input.txt` — original-run bookkeeping, not read by
  the current pipeline.
- `dog10k_mapping_start_at_bams.sh` — original-author reference script (pre-conversion), kept for
  provenance.
