# 3. Merging & filtering

Nextflow workflow that builds the GWAS-ready merged/filtered PLINK datasets. Full detail in
[`code/02_data_processing/3_Merging_filtering/README.md`](../../../code/02_data_processing/3_Merging_filtering/README.md).

## Files

| File | Description |
|---|---|
| [`run_nextflow.sh`](run_nextflow.sh) | Runs the workflow via Nextflow |
| [`params.yml`](params.yml) | Input, frozen-checkpoint, and manual-decision file paths |
| [`nextflow.config`](nextflow.config) | Executor profiles and resource limits |

## Running

Requires [`01_mapping`](../../01_mapping/) and
[`2_Axiom_imputation`](../2_Axiom_imputation/) to have completed first — see `params.yml`.

```bash
export NXF_PROFILE=uppmax
pixi run 00b-mapping-batched
pixi run 02b-axiom-imputation
pixi run 02d-merging-filtering
```
