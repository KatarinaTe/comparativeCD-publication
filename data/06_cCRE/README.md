# data/06_cCRE

Deposited inputs, fetched raw data, and committed intermediate files for
[`code/06_cCRE`](../../code/06_cCRE).

## Live Nextflow pipeline inputs

| File | Consumed by |
|---|---|
| `emissions_state_probabilities.txt`, `emissions_genome_annotations_enrichment.txt` | `EMISSION_STATE_PLOT` (`cCRE.nf`) — deposited-only inputs |
| `dogCD_gws100kb_collapsed.bed`, `SIZE_gwas_clump100kb_regions_for_cCREintersect_clean.txt` | `cre_overlap.nf` — deposited reference copies of `06b-clump-regions-100kb`'s own output (see `cre_overlap/README.md`'s Notes) |

## Live Nextflow pipeline outputs

| File | Produced by |
|---|---|
| `emission_states_plot.{pdf,png}` | `EMISSION_STATE_PLOT` — Fig. 4c |
| `epicdog_tissue_bp_overlaps.txt` | `BUILD_EPICDOG_TISSUE_OVERLAPS` (`cre_overlap.nf`) — bp-level overlap of dogCD/SIZE GWAS regions with pooled brain/other EpicDog tissues, published to `data/results/06_cCRE/cre_overlap/overlaps/`. Its live output (2026-09-23 run) reproduces Supplementary Table 13 exactly; the older copy that sat here (2026-09-17, from pre-final GWAS regions) is in `archive/cCRE_overlap_tables_superseded/` |
| `tissue_specificity.pdf` | `PLOT_TISSUE_SPECIFICITY` — EpicDog odds ratios (Results text, ST13; not a figure panel) |
| `brain_region_bp_overlaps.txt`, `SIZE_brain_region_bp_overlaps.txt` | `BUILD_UU_REGION_BP_OVERLAPS` (`cre_overlap.nf`) — bp-level overlap of dogCD/SIZE GWAS regions with each of the 8 UU brain-region cCRE sets, published to `data/results/06_cCRE/cre_overlap/overlaps/`. Its live output (2026-09-23 run) reproduces the per-region ORs in the Results text (ACG 1.47, cerebellum 1.18); the older copies that sat here are in `archive/cCRE_overlap_tables_superseded/` |
| `region_forest.pdf` | `PLOT_REGION_FOREST` — Fig. 4e |

## Fetched raw data (gitignored — regenerate via the fetch scripts, not committed)

| Directory | Fetched by |
|---|---|
| `chromatin_states_bed/` | `analyses/06_cCRE/fetch-chromatin-states-bed.sh` — 8 per-brain-region UU ChromHMM 9-state BEDs |
| `epicdog_chromatin_states_bed_canfam4/` | `analyses/06_cCRE/fetch-epicdog-chromatin-states-bed.sh` — 11 per-tissue EpicDog CanFam4 13-state BEDs |

## Orphaned — not consumed by anything

| File | Description |
|---|---|
| `brain_region_overlaps.txt` | Per-UU-brain-region **element-count** (not bp) overlap counts — fed a global chi-square/per-region-Fisher/FDR side-analysis in the now-archived `testing_CRE_overlap_260617.R` that isn't part of either published Fig. 4e or 4f (see `cre_overlap/README.md`'s Notes on what was dropped from that script). Kept for now rather than deleted, but nothing reads it. |

## `dogCD_gws100kb_collapsed.bed` also feeds

`compare_UU_epicdog_elements.R` (Fig. 4d) reads this file too, independent of `cre_overlap.nf`.
