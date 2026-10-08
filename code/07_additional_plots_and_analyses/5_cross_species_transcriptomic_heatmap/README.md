# 07. Additional plots and analyses — 5. Cross-species transcriptomic heatmap

Extended Data Fig. 3: heatmap of pairwise Spearman rank correlation coefficients between human
and dog pseudo-bulk expression profiles for the 165 OCD/dogCD network genes, across matched cell
types (neurons vs. glia/other).

**Status: run successfully 2026-09-30 against the real data on the cluster. KatarinaTe has
confirmed these recomputed numbers and this regenerated figure are the ones going forward —
supersede the manuscript's currently-published "Spearman r = 0.80–0.89" figure/claim. See
Verification against the published figure and Notes below.**

**Pipeline engine:** Nextflow (`pixi run 07g-cross-species-transcriptomic-heatmap`)
**Containerized:** Yes — reuses `3_cross_species_expression`'s `cross_species_expression`
container as-is, no new build. It already has `bioconductor-complexheatmap`/`r-circlize`
installed (unused there since 2026-09-22; see
[`../3_cross_species_expression/README.md`](../3_cross_species_expression/README.md)), which is
exactly what this script needs.
**Estimated runtime:** [FILL IN]

---

## Collaborator's original script

[`../dogCD_ED-fig4_human_dog_heatmap.R`](../dogCD_ED-fig4_human_dog_heatmap.R) — dropped in place
2026-09-30, same convention as `3_cross_species_expression`'s
[`network_gene_expression_enrichment_plots_code_updated.r`](../network_gene_expression_enrichment_plots_code_updated.r):
kept loose at the top level of `07_additional_plots_and_analyses/`, content unmodified.
[`bin/cross_species_transcriptomic_heatmap.R`](bin/cross_species_transcriptomic_heatmap.R) is the
literal, CLI-args reproduction — see its own header comment for exactly what changed and why.

## Input Data

| Parameter | File | Source |
|-----------|------|--------|
| `human_dog_rds` | `human_dog_5regions.rds` (4.8 GB) | Downloaded 2026-09-30 via `gdown` from a Google Drive link to `/proj/uppmax2025-2-49/comparativeCD/data/07_additional_plots_and_analyses/human_dog_5regions.rds` on Pelle — not yet re-deposited on Figshare (see Notes). Collaborator's own description: "I generated it from the integrated human-dog-mouse single cell Seurat object I made, after dropping all mouse cells from the dataset." |
| `seed_gene_list` | [`data/05_network/ocd_ccd_systemsmap_genes.txt`](../../../data/05_network/ocd_ccd_systemsmap_genes.txt) (165 genes) | Already deposited — identical to `dogCD_OCD_network`, hardcoded a second time with no visible source in the collaborator's original script; read from this file instead of a third hardcoded copy. |

---

## Cell-type abbreviation mapping (`abbrev_map`)

Referenced in the collaborator's original script but never defined there, and not found anywhere
else in the repo. Resolved 2026-09-30 by reading the real embedded heatmap image directly
(`image2.png` inside `NATURE_Extended_Data_Figures_260929.docx` — the actual published Extended
Data Fig. 3) axis-by-axis, then matching each short axis label back to the real
`supercluster_term` values pulled from `human_dog_5regions.rds` itself
(`unique(dh_5regions@meta.data$supercluster_term)`, 26 values total; 20 are plotted).

The heatmap splits into two blocks — **Neurons** and **Glia/Vascular/Other** — same split the
script uses for `row_split`/`column_split` via the `neurons` vector below.

**Neurons (11)**

| Axis label (`cell_abbrevs`) | `supercluster_term` |
|---|---|
| UL-tel | Upper-layer intratelencephalic |
| DL-tel | Deep-layer intratelencephalic |
| DL-cor-6b | Deep-layer corticothalamic and 6b |
| DL-NP | Deep-layer near-projecting |
| Th-ex | Thalamic excitatory |
| Cer-ex | Upper rhombic lip |
| Cer-i | Cerebellar inhibitory |
| Mid-i | Midbrain-derived inhibitory |
| CGE-inter | CGE interneuron |
| MGE-inter | MGE interneuron |
| LAMP5 | LAMP5-LHX6 and Chandelier |

**Glia / Vascular / Other (9)**

| Axis label (`cell_abbrevs`) | `supercluster_term` |
|---|---|
| Astro | Astrocyte |
| Bergmann | Bergmann glia |
| Ependymal | Ependymal *(not abbreviated)* |
| Fibro | Fibroblast |
| Micro | Microglia |
| Oligo | Oligodendrocyte |
| OP | Oligodendrocyte precursor |
| COP | Committed oligodendrocyte precursor |
| Vascular | Vascular *(not abbreviated)* |

**`Cer-ex` = "Upper rhombic lip"** is the one entry not self-evident from the image alone (the
data has no `supercluster_term` literally called "Cerebellar excitatory"). Resolved by
cross-referencing this repo's own already-converted, manuscript-verified sibling script,
[`3_cross_species_expression/bin/cross_species_expression.R`](../3_cross_species_expression/bin/cross_species_expression.R),
whose own `neurons` vector independently lists `"URL"` as its abbreviation for the same
underlying `supercluster_term` — the developmental lineage that produces excitatory cerebellar
granule neurons, in the Siletti atlas's own terminology. Two independent script contexts landing
on the same underlying cell population is strong corroboration, though the abbreviation string
itself ("Cer-ex" vs. "URL") differs between the two scripts/figures.

**6 of the 26 `supercluster_term` values are not plotted in this heatmap at all** — consistent
with the "20 plotted" KatarinaTe confirmed directly against the published figure. The real
per-species cell counts (from `table(dh_5regions$species, dh_5regions$supercluster_term)`, printed
by the script) show why: each of the 6 has essentially no cells in one of the two species, making a
cross-species correlation for that type meaningless:

| `supercluster_term` | Dog cells | Human cells |
|---|---|---|
| Choroid plexus | 399 | 0 |
| Splatter | 0 | 1,491 |
| Miscellaneous | 0 | 301 |
| Medium spiny neuron | 0 | 87 |
| Eccentric medium spiny neuron | 0 | 68 |
| Amygdala excitatory | 0 | 18 |

---

## Verification against the published figure

Confirmed 2026-09-30: **all 400 cells of the 20×20 human×dog correlation matrix match the
published Extended Data Fig. 3 exactly**, to 2 decimal places — checked cell-by-cell against the
real embedded heatmap image (`image2.png`), not just the diagonal. `abbrev_map`, the
`AverageExpression(..., slot = "data")` call, and the DAOA gene exclusion (below) together
reproduce the published figure exactly; nothing left to reconcile.

One thing worth noting for anyone comparing against the manuscript's Results text specifically:
the text states "Spearman r = 0.80–0.89" for matching neuronal cell types, but 3 of the 11 neuron
diagonal values are marginally outside that stated range — Cer-i (0.80, borderline), Mid-i
(0.79), and Cer-ex (0.78). This isn't a reproduction discrepancy: the published figure's own
printed diagonal shows the exact same three values, so the text's "0.80–0.89" was always a
slightly rounded characterization of the true 0.78–0.89 spread visible in the figure itself, not
a strict per-cell floor.

---

## Proposed figure-legend addendum

For KatarinaTe to fold into the actual Extended Data Fig. 3 legend, in the same style as its
existing abbreviation glossary ("UL = upper layer; DL = deep layer; ..."), which currently only
defines the 11 neuron-block abbreviations:

> Astro = astrocyte; Bergmann = Bergmann glia; COP = committed oligodendrocyte precursor;
> Ependymal = ependymal; Fibro = fibroblast; Micro = microglia; Oligo = oligodendrocyte; OP =
> oligodendrocyte precursor; Vascular = vascular. Cer-ex denotes the upper rhombic lip lineage,
> the developmental source of excitatory cerebellar granule neurons in the Siletti atlas
> nomenclature. Of the 165 OCD/dogCD network genes, DAOA was not detected in this dataset's RNA
> assay and is excluded from the correlations shown (n = 164 genes). Six supercluster types
> present in the source data are not shown, each represented by essentially no cells in one of
> the two species: choroid plexus (399 dog, 0 human), splatter (0 dog, 1,491 human),
> miscellaneous (0 dog, 301 human), medium spiny neuron (0 dog, 87 human), eccentric medium
> spiny neuron (0 dog, 68 human), and amygdala excitatory (0 dog, 18 human).

---

## Notes

- **`col_fun_cor`** (the heatmap's colour ramp) was also referenced but never defined in the
  original script, and — unlike `abbrev_map` — isn't recoverable from the data at all (it's a
  plotting choice, not biological data). Approximated in `bin/` with the standard `viridis`
  palette's own well-known anchor hex colours, in their normal order (low correlation = dark
  purple/blue, high correlation = yellow — matching the published heatmap) rather than pulling in
  a new package dependency for one colour ramp. Briefly reversed 2026-09-30, switched back
  2026-10-01 per KatarinaTe once compared against the original figure. Doesn't affect any of the
  actual correlation numbers, only the heatmap's rendered colours — revisit only if exact
  colour-for-colour figure reproduction matters later.
- **Data availability**: `human_dog_5regions.rds` is fetched from Figshare, DOI
  `10.6084/m9.figshare.33787129` (same article as `Siletti_human_subset_normalised.rds`, used by
  the sibling `3_cross_species_expression` stage — this file was added to that same deposit
  rather than getting its own). First checked 2026-09-30, when this DOI only listed the Siletti
  file yet; confirmed live 2026-10-05, checked directly via Figshare's API
  (`api.figshare.com/v2/articles/33787129`), now listing both files.
  `07f-fetch-human-dog-5regions` (`analyses/.../fetch-human-dog-5regions.sh`) fetches it, wired
  into `07g-cross-species-transcriptomic-heatmap`'s `depends-on`. Confirmed byte-identical
  (same size, 5,115,368,146 bytes, and md5, `3812a1c0743bd506359fb77f1e158391`) to the
  cluster-local copy fetched manually via `gdown` on 2026-09-30, which is what the verification
  below actually ran against — so that verification and the now-fetched copy are provably the
  same data, not just presumed to be.
- Nextflow-converted 2026-09-30, following the same pattern as `3_cross_species_expression` (see
  [`archive/conversion_notes/07_additional_plots_and_analyses.md`](../../../archive/conversion_notes/07_additional_plots_and_analyses.md)
  and its sibling stages' notes for the general conversion approach used throughout this repo).
  Reuses that stage's exact container — no new build needed.

---

## How to run

```bash
export NXF_PROFILE=uppmax
pixi run 07g-cross-species-transcriptomic-heatmap
```

Runs `bin/cross_species_transcriptomic_heatmap.R` inside the `cross_species_expression` container
via Nextflow, reading paths from
[`analyses/07_additional_plots_and_analyses/5_cross_species_transcriptomic_heatmap/params.yml`](../../../analyses/07_additional_plots_and_analyses/5_cross_species_transcriptomic_heatmap/params.yml)
(currently the cluster-local `human_dog_rds` path — see Notes). Outputs publish to
`data/results/07_additional_plots_and_analyses/5_cross_species_transcriptomic_heatmap/`.

Stub test:

```bash
NXF_PROFILE=stub_test ./run_nextflow.sh   # from the analyses/ directory above
```

### Standalone (outside Nextflow), for ad hoc runs

Same container, run directly via Apptainer/Singularity (substitute `singularity` for `apptainer`
below if that's what's available instead):

```bash
apptainer exec /gorilla/proj/uppmax2025-2-49/comparativeCD/env/apptainer_cache/cross_species_expression.sif \
    Rscript /gorilla/proj/uppmax2025-2-49/comparativeCD/code/07_additional_plots_and_analyses/5_cross_species_transcriptomic_heatmap/bin/cross_species_transcriptomic_heatmap.R \
    /gorilla/proj/uppmax2025-2-49/comparativeCD/data/07_additional_plots_and_analyses/human_dog_5regions.rds \
    /gorilla/proj/uppmax2025-2-49/comparativeCD/data/05_network/ocd_ccd_systemsmap_genes.txt \
    <output directory>
```

(Absolute, real — not symlinked — paths throughout: `/proj` on this cluster is a symlink to
`/gorilla/proj`, which broke both plain relative paths and `--bind /proj:/proj` when this was
first run standalone. See git history for the debugging trail if this resurfaces.)

---

## Process Map

| # | Process | Container | Description |
|---|---------|-----------|-------------|
| 1 | `CROSS_SPECIES_TRANSCRIPTOMIC_HEATMAP` | `cross_species_expression` | Runs `bin/cross_species_transcriptomic_heatmap.R` end to end — average expression, correlation matrix, heatmap PDF |

---

## Output Data

| Path | Description |
|------|-------------|
| `human_dog_correlation_by_celltype.csv` | full-precision Spearman r + p-value, one row per cell type with data in both species (superset of what's plotted) |
| `human_dog_correlation_matrix.rds` | cached Spearman correlation matrix (human x dog, the 20 matched/ordered cell types in the published figure) |
| `human_dog_heatmap.pdf` | Extended Data Fig. 3 |
