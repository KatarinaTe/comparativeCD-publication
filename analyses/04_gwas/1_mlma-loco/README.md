# 1. MLMA-LOCO (factors + SIZE)

Nextflow workflow for Fig. 1c's mlma-loco GWAS (F1/F2/F3 + SIZE). Full detail in
[`code/04_gwas/1_mlma-loco/README.md`](../../../code/04_gwas/1_mlma-loco/README.md).

## Files

| File | Description |
|---|---|
| [`run_nextflow.sh`](run_nextflow.sh) | Runs the workflow via Nextflow |
| [`params.yml`](params.yml) | QC6 plink/phenotype directory paths, gene-range file |
| [`nextflow.config`](nextflow.config) | Executor profiles and resource limits |

## Running

Requires [`3_Merging_filtering`](../../02_data_processing/3_Merging_filtering/) to have
completed first — see `params.yml`.

```bash
export NXF_PROFILE=uppmax
pixi run 02d-merging-filtering
pixi run 04a-mlma-loco
```
