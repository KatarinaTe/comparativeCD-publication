# 4. QC signals (IGV)

Nextflow workflow that preps inputs for the manual IGV signal-QC step. **The IGV inspection itself
is manual, not run by this workflow.** Full detail in
[`code/04_gwas/4_qc_signals_IGV/README.md`](../../../code/04_gwas/4_qc_signals_IGV/README.md).

## Files

| File | Description |
|---|---|
| [`run_nextflow.sh`](run_nextflow.sh) | Runs the workflow via Nextflow |
| [`params.yml`](params.yml) | Fine-mapping results dir, BAM dir, QC5 plink dir, crosswalk fam, individual-selection tuning |
| [`nextflow.config`](nextflow.config) | Executor profiles and resource limits |

## Running

Requires [`3_finemap_susie`](../3_finemap_susie/) and
[`3_Merging_filtering`](../../02_data_processing/3_Merging_filtering/) to have completed first,
and LowPass BAMs materialized (`pixi run 00b-mapping-batched`) — see `params.yml`.

```bash
export NXF_PROFILE=uppmax
pixi run 04c-finemap-susie
pixi run 04d-qc-signals-igv
```

Re-run with a different `--round` value (in `params.yml`) for a fresh individual-selection batch.
