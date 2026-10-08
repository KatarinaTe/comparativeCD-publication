# 5. Plotting

Nextflow workflow for Fig. 1c's Manhattan + QQ plots. Full detail in
[`code/04_gwas/5_plotting/README.md`](../../../code/04_gwas/5_plotting/README.md).

## Files

| File | Description |
|---|---|
| [`run_nextflow.sh`](run_nextflow.sh) | Runs the workflow via Nextflow |
| [`params.yml`](params.yml) | `1_mlma-loco`/`2_polmm` published `gwas/` directory paths |
| [`nextflow.config`](nextflow.config) | Executor profiles and resource limits |

## Running

Requires [`1_mlma-loco`](../1_mlma-loco/) and [`2_polmm`](../2_polmm/) to have completed first —
see `params.yml`.

```bash
export NXF_PROFILE=uppmax
pixi run 04a-mlma-loco
pixi run 04b-polmm
pixi run 04e-gwas-plotting
```
