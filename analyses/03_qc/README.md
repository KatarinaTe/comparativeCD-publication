# 03. QC

Nextflow workflow for the genotype-concordance QC studies. Full detail in
[`code/03_qc/README.md`](../../code/03_qc/README.md).

## Files

| File | Description |
|---|---|
| [`run_nextflow.sh`](run_nextflow.sh) | Runs the workflow via Nextflow with `params.yml` |
| [`params.yml`](params.yml) | Reference file paths, cohort input paths, and sample/pair samplesheet paths |
| [`nextflow.config`](nextflow.config) | Executor profiles and resource limits |

## Running

Requires [`mapping`](../01_mapping/) to have completed first (studies 1 and 3 need its HighPass
outputs; studies 1 and 2 need its LowPass outputs). Study 2 additionally needs
`02_data_processing/2_Axiom_imputation` to have completed (its published `vcf/` output feeds
`axiom_lowpass_vcf_dir` in `params.yml`).

```bash
export NXF_PROFILE=uppmax
pixi run 00b-mapping-batched
pixi run 03a-qc
```
