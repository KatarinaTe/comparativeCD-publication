# 5. Cross-species transcriptomic heatmap

Nextflow workflow for Extended Data Fig. 3's cross-species transcriptomic concordance heatmap.
Verified 2026-09-30 against the published figure: all 400 cells of the recomputed 20×20 human×dog
correlation matrix match exactly. Full detail in
[`code/07_additional_plots_and_analyses/5_cross_species_transcriptomic_heatmap/README.md`](../../../code/07_additional_plots_and_analyses/5_cross_species_transcriptomic_heatmap/README.md).

**`human_dog_rds` is fetched from Figshare** (DOI
[10.6084/m9.figshare.33787129](https://doi.org/10.6084/m9.figshare.33787129), activated
2026-10-05 once the deposit went live) by `07f-fetch-human-dog-5regions`, now wired into
`07g-cross-species-transcriptomic-heatmap`'s `depends-on`. Confirmed byte-identical to the
cluster-local copy (fetched manually via `gdown` on 2026-09-30) the 2026-09-30 verification run
actually used — same size and md5 — so that verification and this fetched copy are provably the
same data.

## Files

| File | Description |
|---|---|
| [`run_nextflow.sh`](run_nextflow.sh) | Runs the workflow via Nextflow |
| [`params.yml`](params.yml) | Path to `human_dog_5regions.rds` (fetched by `07f`, see above) and the deposited 165-gene seed list |
| [`nextflow.config`](nextflow.config) | Executor profiles and resource limits |
| [`fetch-human-dog-5regions.sh`](fetch-human-dog-5regions.sh) | Fetches `human_dog_5regions.rds` from Figshare (10.6084/m9.figshare.33787129) |

## Running

```bash
export NXF_PROFILE=uppmax
pixi run 07g-cross-species-transcriptomic-heatmap
```
