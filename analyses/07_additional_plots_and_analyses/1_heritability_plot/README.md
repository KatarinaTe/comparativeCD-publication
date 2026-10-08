# 1. Heritability plot

Nextflow workflow for Fig. 1b's heritability boxplot. No pipeline prerequisites. Full detail in
[`code/07_additional_plots_and_analyses/1_heritability_plot/README.md`](../../../code/07_additional_plots_and_analyses/1_heritability_plot/README.md).

## Files

| File | Description |
|---|---|
| [`run_nextflow.sh`](run_nextflow.sh) | Runs the workflow via Nextflow |
| [`params.yml`](params.yml) | Path to `NATURE_SupplementaryTables_261006.xlsx` |
| [`nextflow.config`](nextflow.config) | Executor profiles and resource limits |

## Running

```bash
export NXF_PROFILE=uppmax
pixi run 07a-heritability-plot
```
