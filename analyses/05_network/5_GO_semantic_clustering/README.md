# 5. GO semantic clustering

Nextflow workflow for the GO semantic-clustering → fold-enrichment → community-tileplot chain
(Fig. 2c-e, 3b, 3d-e). Full detail in
[`code/05_network/5_GO_semantic_clustering/README.md`](../../../code/05_network/5_GO_semantic_clustering/README.md).

## Files

| File | Description |
|---|---|
| [`run_nextflow.sh`](run_nextflow.sh) | Runs the workflow via Nextflow |
| [`params.yml`](params.yml) | Deposited per-disease GO-enrichment TSV paths, cached GPT labels, GO annotation file |
| [`nextflow.config`](nextflow.config) | Executor profiles and resource limits |

## Running

No pipeline prerequisites — reads deposited GO-enrichment TSVs directly (does not depend on
`1_NetColoc`/`2_CrossSpeciesBMI`'s own outputs).

```bash
export NXF_PROFILE=uppmax
pixi run 05c-go-semantic-clustering
```
