# 06. cCRE — Clump regions (±100kb)

Nextflow workflow that builds ±100kb GWAS clump regions from `1_mlma-loco`'s and `2_polmm`'s
output. Full detail in
[`code/06_cCRE/clump_regions_100kb/README.md`](../../../code/06_cCRE/clump_regions_100kb/README.md).

## Files

| File | Description |
|---|---|
| [`run_nextflow.sh`](run_nextflow.sh) | Runs the workflow via Nextflow |
| [`params.yml`](params.yml) | `1_mlma-loco`'s and `2_polmm`'s published `clumping/`/`gwas/` directory paths |
| [`nextflow.config`](nextflow.config) | Executor profiles and resource limits |

## Running

Requires [`1_mlma-loco`](../../04_gwas/1_mlma-loco/) and [`2_polmm`](../../04_gwas/2_polmm/) to
have completed first — see `params.yml`.

```bash
export NXF_PROFILE=uppmax
pixi run 04a-mlma-loco
pixi run 04b-polmm
pixi run 06b-clump-regions-100kb
```
