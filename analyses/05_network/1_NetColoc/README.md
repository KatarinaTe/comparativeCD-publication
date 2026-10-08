# 1. NetColoc

Nextflow workflow for NetColoc network-propagation z-scores over PCNet2.0. Full detail in
[`code/05_network/1_NetColoc/README.md`](../../../code/05_network/1_NetColoc/README.md).

## Files

| File | Description |
|---|---|
| [`run_nextflow.sh`](run_nextflow.sh) | Runs the workflow via Nextflow |
| [`params.yml`](params.yml) | NDEx interactome UUID, seed-gene directory |
| [`nextflow.config`](nextflow.config) | Executor profiles and resource limits |

## Running

```bash
export NXF_PROFILE=uppmax
pixi run 05a-netcoloc
```

No pipeline prerequisites — reads deposited seed-gene files and fetches PCNet2.0 from NDEx
directly.
