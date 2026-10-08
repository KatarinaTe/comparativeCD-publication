# 2. Axiom imputation

Nextflow workflow for Axiom-to-Dog10K genotype imputation. Full detail in
[`code/02_data_processing/2_Axiom_imputation/README.md`](../../../code/02_data_processing/2_Axiom_imputation/README.md).

## Files

| File | Description |
|---|---|
| [`run_nextflow.sh`](run_nextflow.sh) | Runs the workflow via Nextflow |
| [`params.yml`](params.yml) | Input, reference, and tool jar paths |
| [`nextflow.config`](nextflow.config) | Executor profiles and resource limits |

## Running

Requires [`fetch-raw-data`](../../00_fetch-raw-data/) and `fetch-axiom-jars` to have completed
first, and `Affy_merged` (the raw Axiom genotypes) to be supplied — see `params.yml`.

```bash
export NXF_PROFILE=uppmax
pixi run 00a-fetch-reference-data
pixi run 02a-fetch-axiom-jars
pixi run 02b-axiom-imputation
```
