# 4. Plotting

Nextflow workflow for `GO_TERM_COMPARISON` only (Fig. comparing human-only vs. cross-species GO
enrichment) — **this directory also holds `Fig2_network-overlap.R` (Fig. 2b), which is a separate,
standalone script not wired into this Nextflow pipeline**; see
[`code/05_network/4_Plotting/README.md`](../../../code/05_network/4_Plotting/README.md) for why.

## Files

| File | Description |
|---|---|
| [`run_nextflow.sh`](run_nextflow.sh) | Runs the Nextflow workflow (`GO_TERM_COMPARISON` only) |
| [`params.yml`](params.yml) | Deposited GO-enrichment TSV directory |
| [`nextflow.config`](nextflow.config) | Executor profiles and resource limits |

## Running

No pipeline prerequisites — reads deposited GO-enrichment TSVs directly, does not depend on
`1_NetColoc`/`2_CrossSpeciesBMI`'s own outputs.

```bash
export NXF_PROFILE=uppmax
pixi run 05c-go-term-comparison
```

`Fig2_network-overlap.R` is run manually (`Rscript Fig2_network-overlap.R`), not via `pixi run`.
