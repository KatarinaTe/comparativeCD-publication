# 1. Survey data

Nextflow workflow for the dogCD survey's EFA/IRT factor analysis. Full detail in
[`code/02_data_processing/1_Survey_data/README.md`](../../../code/02_data_processing/1_Survey_data/README.md).

## Files

| File | Description |
|---|---|
| [`run_nextflow.sh`](run_nextflow.sh) | Runs the workflow via Nextflow |
| [`params.yml`](params.yml) | Response table and survey metadata paths |
| [`nextflow.config`](nextflow.config) | Executor profiles and resource limits |

## Running

```bash
export NXF_PROFILE=uppmax
pixi run 02c-survey-data
```
