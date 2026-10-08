# 2. GWAS catalog overlap

Nextflow workflow reproducing the "12 of 27" GWAS-catalog claim's mechanical join. No pipeline
prerequisites — reads frozen inputs. Full detail in
[`code/07_additional_plots_and_analyses/2_gwas_catalog_overlap/README.md`](../../../code/07_additional_plots_and_analyses/2_gwas_catalog_overlap/README.md).

## Files

| File | Description |
|---|---|
| [`run_nextflow.sh`](run_nextflow.sh) | Runs the workflow via Nextflow |
| [`params.yml`](params.yml) | Paths to the archived GWAS-catalog originals + the frozen gene-list stand-in |
| [`nextflow.config`](nextflow.config) | Executor profiles and resource limits |

## Running

```bash
export NXF_PROFILE=uppmax
pixi run 07b-gwas-catalog-overlap
```

The `rating`/community-membership/`comment` columns of the frozen gene list remain manually
curated — see `REPRODUCIBILITY_AUDIT.md` for why that part stays out of scope for automation.
