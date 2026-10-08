# 4. Power analysis (Fig 1d)

Nextflow workflow reproducing the dogCD-vs-human-OCD GWAS power comparison. No pipeline
prerequisites — reads frozen inputs. Full detail in
[`code/07_additional_plots_and_analyses/4_power_analysis/README.md`](../../../code/07_additional_plots_and_analyses/4_power_analysis/README.md).

## Files

| File | Description |
|---|---|
| [`run_nextflow.sh`](run_nextflow.sh) | Runs the workflow via Nextflow |
| [`params.yml`](params.yml) | Paths to the frozen dogCD/human-OCD GWAS input tables |
| [`nextflow.config`](nextflow.config) | Executor profiles and resource limits |

## Running

```bash
export NXF_PROFILE=uppmax
pixi run 07e-power-analysis
```
