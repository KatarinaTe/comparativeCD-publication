# 06. cCRE — emission-state-plot

Nextflow workflow for Fig. 4c's emission-probability heatmap. No pipeline prerequisites — reads
deposited-only inputs. Full detail in [`code/06_cCRE/README.md`](../../code/06_cCRE/README.md).

## Files

| File | Description |
|---|---|
| [`run_nextflow.sh`](run_nextflow.sh) | Runs the workflow via Nextflow |
| [`params.yml`](params.yml) | Deposited emission-probability/genome-annotation-enrichment file paths |
| [`nextflow.config`](nextflow.config) | Executor profiles and resource limits |
| [`fetch-chromatin-states-bed.sh`](fetch-chromatin-states-bed.sh) | Fetches the 8 per-brain-region ChromHMM 9-state BEDs (Figshare, checksum-verified) |
| [`fetch-epicdog-chromatin-states-bed.sh`](fetch-epicdog-chromatin-states-bed.sh) | Fetches the 11 per-tissue EpicDog CanFam4 13-state BEDs (Figshare, checksum-verified) |

## Running

```bash
export NXF_PROFILE=uppmax
pixi run 06a-emission-state-plot
```

The two `fetch-*.sh` scripts are run manually, not via `pixi run` — see
[`code/06_cCRE/README.md`](../../code/06_cCRE/README.md) for what consumes their output: the
loose Fig. 4d comparison script, and [`cre_overlap/`](cre_overlap) (Fig. 4e, its own Nextflow
substage, `pixi run 06c-cre-overlap`).
