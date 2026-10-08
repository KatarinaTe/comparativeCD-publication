# 3. Cross-species expression

Nextflow workflow for Fig. 4a-b's gene-set enrichment vs. the Siletti snRNA-seq atlas. Needs the
deposited Siletti subset fetched first — wired in via `07d-cross-species-expression`'s
`depends-on = ["07c-fetch-siletti-subset"]`, so `pixi run 07d-cross-species-expression` alone
fetches it automatically if not already present. Full detail in
[`code/07_additional_plots_and_analyses/3_cross_species_expression/README.md`](../../../code/07_additional_plots_and_analyses/3_cross_species_expression/README.md).

## Files

| File | Description |
|---|---|
| [`run_nextflow.sh`](run_nextflow.sh) | Runs the workflow via Nextflow |
| [`params.yml`](params.yml) | Path to the fetched Siletti subset `.rds` |
| [`nextflow.config`](nextflow.config) | Executor profiles and resource limits |
| [`fetch-siletti-subset.sh`](fetch-siletti-subset.sh) | Fetches `Siletti_human_subset_normalised.rds` (Figshare, 2.58 GB, md5-verified) — `pixi run 07c-fetch-siletti-subset` |

## Running

```bash
export NXF_PROFILE=uppmax
pixi run 07d-cross-species-expression   # fetches the Siletti subset first (07c) if needed, then runs
```

The script that produced this `.rds` subset from the full Siletti et al. 2023 atlas cannot be
provided (per the manuscript author) — reproducibility stops at this deposited file. See
`REVIEW_AND_SUGGESTIONS.md`'s Fig 4a-b entry.
