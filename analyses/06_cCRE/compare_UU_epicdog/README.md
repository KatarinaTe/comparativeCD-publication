# 06. cCRE — Compare UU vs. EpicDog elements

Nextflow workflow that compares element counts between UU's own cCRE atlas and EpicDog's, for
Fig. 4d. Full detail in
[`code/06_cCRE/compare_UU_epicdog/README.md`](../../../code/06_cCRE/compare_UU_epicdog/README.md).

## Files

| File | Description |
|---|---|
| [`run_nextflow.sh`](run_nextflow.sh) | Runs the workflow via Nextflow |
| [`params.yml`](params.yml) | EpicDog + UU fetched-BED paths (8 UU per-region BEDs, 2 EpicDog BEDs) |
| [`nextflow.config`](nextflow.config) | Executor profiles and resource limits |

## Running

Requires both sides fetched first — see `params.yml`.

```bash
export NXF_PROFILE=uppmax
pixi run 06d-fetch-epicdog-bed
pixi run 06e-fetch-chromatin-states-bed
pixi run 06f-compare-uu-epicdog
```
