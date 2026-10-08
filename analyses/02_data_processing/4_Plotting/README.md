# 4. Plotting

Nextflow workflow for terminal diagnostic plots over `3_Merging_filtering`'s outputs. Full detail
in [`code/02_data_processing/4_Plotting/README.md`](../../../code/02_data_processing/4_Plotting/README.md).

## Files

| File | Description |
|---|---|
| [`run_nextflow.sh`](run_nextflow.sh) | Runs the workflow via Nextflow |
| [`params.yml`](params.yml) | Input file/directory paths |
| [`nextflow.config`](nextflow.config) | Executor profiles and resource limits |

## Running

Requires [`3_Merging_filtering`](../3_Merging_filtering/) to have completed first — see
`params.yml`.

```bash
export NXF_PROFILE=uppmax
pixi run 02d-merging-filtering
pixi run 02e-plotting
```
