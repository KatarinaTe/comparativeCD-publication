# 2. CrossSpeciesBMI

Nextflow workflow for the dogCD/OCD conserved network and colocalization tests (Fig. 1e-i, 2a, 3a).
Full detail in
[`code/05_network/2_CrossSpeciesBMI/README.md`](../../../code/05_network/2_CrossSpeciesBMI/README.md).

## Files

| File | Description |
|---|---|
| [`run_nextflow.sh`](run_nextflow.sh) | Runs the workflow via Nextflow |
| [`params.yml`](params.yml) | `1_NetColoc` output paths + deposited seed-gene/control/hierarchy/MGD file paths |
| [`nextflow.config`](nextflow.config) | Executor profiles and resource limits |

## Running

Requires [`1_NetColoc`](../1_NetColoc/) to have completed first — see `params.yml`.

```bash
export NXF_PROFILE=uppmax
pixi run 05a-netcoloc
pixi run 05b-cross-species-bmi
```
