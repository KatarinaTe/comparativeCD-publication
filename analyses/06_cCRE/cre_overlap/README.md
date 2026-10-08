# 06. cCRE — cCRE overlap (Fig. 4e)

Nextflow workflow for Fig. 4e (UU per-region) and the EpicDog tissue-specificity odds ratios
(Results text, Supplementary Table 13), sharing
GWAS-region inputs. Full detail in
[`code/06_cCRE/cre_overlap/README.md`](../../../code/06_cCRE/cre_overlap/README.md).

## Files

| File | Description |
|---|---|
| [`run_nextflow.sh`](run_nextflow.sh) | Runs the workflow via Nextflow |
| [`params.yml`](params.yml) | Chromatin-state directory and GWAS-region file paths |
| [`nextflow.config`](nextflow.config) | Executor profiles and resource limits |

## Running

Requires the EpicDog chromatin-state BEDs fetched first
(`analyses/06_cCRE/fetch-epicdog-chromatin-states-bed.sh`) and `06b-clump-regions-100kb` run first
(pulled in automatically via `depends-on`) — see `params.yml`.

```bash
export NXF_PROFILE=uppmax
pixi run 06d-fetch-epicdog-bed
pixi run 06c-cre-overlap
```
