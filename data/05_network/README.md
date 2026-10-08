# data/05_network

Deposited/frozen inputs for [`code/05_network`](../../code/05_network). Mostly loose files at
this top level (no clean per-substage split, since these predate the Nextflow conversion) plus a
few subdirectories.

## Seed gene lists (`1_NetColoc` input)

`OCDgenes_v4.txt`, `CCDgenes_Oct25.txt`, `human_DEP_noMHC_seed_genes_header.txt`,
`human_schiz_seed_genes_header.txt` — the 4 gene sets NetColoc z-scores. `*_noheader.txt` variants
and the non-`_header` DEP/SCH files are alternate/older versions not read by the current pipeline.
Control-trait seed genes (`DogHeight_seed_genes.txt`, `human_height2_seed_genes.txt`,
`human_RA_seed_genes.txt`, `human_type2diabetes_seed_genes.txt`) feed the control-comparison
analysis described in `2_CrossSpeciesBMI`.

## Deposited GO/MGD enrichment results (`2_CrossSpeciesBMI` / `5_GO_semantic_clustering` input)

Per-disease-pairing files, one set per pairing (`OCD_ONLY`, `CCD_ONLY`/dogCD, `DEP_ONLY`,
`SCH_ONLY`, plus the 3 cross-species pairings `Compulsive_hierachy_full_GO_enrichment_251023.tsv`
[OCD/dogCD], `DEP_CCD_hierachy_full_GO_enrichment.tsv`, `SCH_CCD_hierachy_full_GO_enrichment.tsv`):
`*_hierachy_full_GO_enrichment.tsv` (GO term enrichment) and
`*_hierarchy_full_MGD_enrichment_results.tsv` (mouse-phenotype enrichment).
`fixedOCD_ONLY_hierarchy_full_MGD_enrichment_results.csv` is `reformat_ocd_only_mgd_results.py`'s
normalized output — see `2_CrossSpeciesBMI/README.md` for why `OCD_ONLY`'s deposit needed
reformatting.

**Neither the GO nor the MGD enrichment itself is re-executed live by this repo** — the tests are
deterministic, but both depend on external reference databases (g:Profiler's live GO:BP backend;
the MGI/MPO phenotype ontology from `informatics.jax.org`) that aren't pinned/versioned anywhere
here, so a rerun today would query whatever state those databases are in now, not the original
snapshot. See `2_CrossSpeciesBMI/README.md`'s Notes section for detail.

## NetColoc / systems-map outputs (deposited, also reproducible via `1_NetColoc`/`2_CrossSpeciesBMI`)

`z_D1_*_z-scores_pcnet2*.csv` / `z_D2_CCD_Oct25_z-scores_pcnet2_251022.csv` (single-species
z-scores), `*_zcomb_z12*.txt` (cross-species combined z-scores), `*_systems_map.txt` /
`CompulsiveNetwork_hierarchy_data_251022.tsv` / `OCD_ONLY_Network_hierarchy_data.tsv` /
`SCH_CCD_Network_hierarchy_data.tsv` (systems-map hierarchies), `*_systemsmap_genes.txt`
(per-pairing gene lists), `toga_cf4_human_genes.csv` (TOGA dog-human orthologue mapping).

**Two exceptions here — deposited but *not* re-executed live by default**, even though the code to
reproduce both exists (source notebook cells):
- `controlanalyses_260521.csv` (control-trait comparison, Fig. 1i and Supplementary Table 10) —
  see `2_CrossSpeciesBMI/README.md`'s Notes section. Its underlying permutation null is unseeded
  in the original code, so reruns aren't numerically identical to this frozen file;
  `CONTROL_ANALYSIS_PLOT` reloads it rather than re-deriving it.
- `netcoloc_enrichment_df_CCD_OCD_251022.csv` (NPS-threshold sensitivity analysis, Supplementary
  Table 7) — see `1_NetColoc/README.md`'s Notes section. The source cell that produces it is
  explicitly marked "do not rerun" and was never ported into a pipeline process.

## GO-term and gene-basis tables (Supplementary Tables 8 & 11)

`all_networks_GO_enrichment_terms.csv/.xlsx` (all significant GO terms across all 7 networks, one
column identifying which network each term belongs to), `systemsmap_community_GO_enrichment_terms.csv/.xlsx`
(same, restricted to each network's top-level systems-map community, not subcommunities),
`OCD_dogCD_all_subcommunities_GO_enrichment_terms.csv/.xlsx` (all significant GO terms for every
subcommunity — C185/C186/C197/C187 — of the OCD/dogCD network specifically; subcommunities were not
used in semantic clustering for any other network, only the OCD/dogCD top level), and
`gene_level_network_basis_table.csv/.xlsx` (one row per gene, one column per seed/network/hierarchy/
PCNet2.0-membership criterion, plus the same four OCD/dogCD subcommunities) — deposited 2026-09-28,
the source data for Supplementary Tables 8 and 11.

**No generating script is currently committed for these four tables** — unlike ST5/ST6, which have
a standalone script ([`code/04_gwas/3_finemap_susie/bin/build_st5_st6_region_ids.py`](../../code/04_gwas/3_finemap_susie/bin/build_st5_st6_region_ids.py)),
these four were built ad hoc against the same deposited GO-enrichment and network files listed above
and checked into the repo as data-only. Revisit whether a standalone script equivalent to
`build_st5_st6_region_ids.py` is worth committing here too.

**`venn_region_{dogCD,OCD,OCD_dogCD}` columns added 2026-10-05** to
`gene_level_network_basis_table.csv/.xlsx` — a human-readable label of which Fig. 2b Venn region
each gene falls into per comparison (e.g. "Depression + Schizophrenia", "OCD only"), via
[`code/05_network/4_Plotting/bin/add_venn_region_labels.R`](../../code/05_network/4_Plotting/bin/add_venn_region_labels.R)
(`ADD_VENN_REGION_LABELS` process, wired into `05d-network-overlap`). Unlike the four tables above,
this one *does* have a generating script, since it's a derived column (not raw network-basis data)
computed from the exact same `in_hierarchy_*` columns and region logic
`plot_venn_hierarchy.R` already uses to draw the Venns — verified to reproduce the same 7 region
counts per comparison as the published figure (see that script's own header). Pre-existing columns
are untouched byte-for-byte; only the 3 new columns were appended.

## Removed (2026-09-18): orphaned pre-conversion overlap/Venn exploration

14 files (`overlap_dep_ocd*.txt`, `overlap_sch_ocd*.txt`, `overlap_dep_ocdCCD*.txt`,
`overlap_dep_sch*.txt`, `overlap_sch_ocdCCD*.txt`, `venn_dep_ocd_sch.pdf`,
`venn_dep_ocdCCD_sch.pdf`, `overlap_stacked*.pdf`, `overlap_proportions_stacked.pdf`) were checked
and removed. Before deleting: confirmed zero references anywhere in `code/`/`analyses/` *and* the
manuscript text; verified all 8 `.txt` files' gene counts against real
`intersect()`/`setdiff()` computations over the systems-map gene lists still in this repo — 7 of 8
matched exactly (real, correct, just unused leftovers from an earlier manual exploration), and the
8th (`overlap_dep_sch_ocd.txt`) turned out to be a byte-identical duplicate of
`overlap_dep_sch.txt` saved under a misleading name (its filename implied a 3-way overlap; its
content was the plain pairwise one). None were an input to, or evidence used for, anything in the
current pipeline or manuscript.

## Subdirectories

| Directory | Description |
|---|---|
| [`5_GO_semantic_clustering/`](5_GO_semantic_clustering) | Frozen inputs specific to that substage (cached GPT cluster labels, GO annotation file) — see its own contents / `code/05_network/5_GO_semantic_clustering/README.md` |
| [`cCRE/`](cCRE) | Brain-region cCRE overlap counts — superseded by `data/06_cCRE/` after cCRE was split into its own stage; kept here for provenance |
| [`cytoscape/`](cytoscape) | `.cx` network exports from NDEx, the reproducible inputs to the manual `3_Cytoscape` step — see that substage's README |
| [`Figures/`](Figures) | Reference figure outputs |
| [`plotting/`](plotting) | Additional plotting inputs/outputs |
