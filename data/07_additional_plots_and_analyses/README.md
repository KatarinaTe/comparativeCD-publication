# data/07_additional_plots_and_analyses

Deposited/frozen inputs for [`code/07_additional_plots_and_analyses`](../../code/07_additional_plots_and_analyses).

## Live pipeline inputs

| File | Consumed by |
|---|---|
| `gene_list_frozen.csv` | `2_gwas_catalog_overlap` — frozen stand-in for the manually-curated `gene_list` (see that substage's own README and `REPRODUCIBILITY_AUDIT.md`) |
| `gwas_catalog_hits.csv` | `2_gwas_catalog_overlap` — raw GWAS Catalog query results. Moved from `archive/gwas_catalog_enrichment_originals/` (2026-09-21) — it's a real pipeline input, not provenance-only material |
| `name_conversion_v2.tsv` | `2_gwas_catalog_overlap` — curated classification of the 5,342 GWAS Catalog reported traits (short name, category, include flag), made by co-authors of this study; local copy of the former "name conversion" sheet. Moved from `archive/` alongside the above |
| `gene_list_gwas_catalog_overlap_with_region.csv` | `2_gwas_catalog_overlap` — one manual merged-result export, used for the crosscheck output. Moved from `archive/` alongside the above |
| `Siletti_human_subset_normalised.rds` | `3_cross_species_expression` — fetched by `fetch-siletti-subset.sh`, not committed (2.58 GB) |
| `human_cell_specificity_scores.rds` | `3_cross_species_expression` — derived by that substage itself from the Siletti subset |
| `ST5.dogCD.top_gwas.tsv` | `4_power_analysis` — frozen Supplementary Table 5 content (per-region clumped top hits, all 17 dogCD traits). Moved from `data/Elinor_simple_power_analysis/` (2026-09-24) — the rest of that bundle (original script, pre-computed outputs, figure-legend docx) is in `archive/elinor_power_analysis_originals/` |
| `Strom.Table1.tsv` | `4_power_analysis` — human OCD GWAS genome-wide-significant loci (Strom et al. 2025 Table 1). Moved alongside the above |

## Removed (2026-09-18): `elinors_analysis/`

Contained `.RData`/`.Rhistory`/`.Rproj.user` and several `cache_go_clusters*.rds` files —
leftover interactive-R-session cache from the collaborator's original working directory, distinct
from (and not to be confused with) the real, referenced
`data/05_network/5_GO_semantic_clustering/elinors_analysis/` cache directory that
`5_GO_semantic_clustering`'s `nextflow.config` actually points at. Not referenced by any current
script or the manuscript text; also not git-tracked at all (`.RData` is gitignored, the rest was
apparently never committed) — removed directly from disk, no git history to preserve.
