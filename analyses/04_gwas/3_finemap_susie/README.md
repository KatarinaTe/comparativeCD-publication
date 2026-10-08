# 3. Fine-mapping (SuSiE via PolyFun)

Nextflow workflow for PolyFun/SuSiE fine-mapping across all 18 GWAS'd phenotypes. Full detail in
[`code/04_gwas/3_finemap_susie/README.md`](../../../code/04_gwas/3_finemap_susie/README.md).

## Files

| File | Description |
|---|---|
| [`run_nextflow.sh`](run_nextflow.sh) | Runs the workflow via Nextflow |
| [`params.yml`](params.yml) | QC5 plink dir, both GWAS stages' sumstats dirs, deposited SIZE/STUCK/ALLFAM `.fam` paths |
| [`nextflow.config`](nextflow.config) | Executor profiles and resource limits |

## Running

Requires [`3_Merging_filtering`](../../02_data_processing/3_Merging_filtering/),
[`1_mlma-loco`](../1_mlma-loco/), and [`2_polmm`](../2_polmm/) to have completed first — see
`params.yml`.

```bash
export NXF_PROFILE=uppmax
pixi run 02d-merging-filtering
pixi run 04a-mlma-loco
pixi run 04b-polmm
pixi run 04c-finemap-susie
```
