# 06. cCRE — Clump regions (±100kb)

Builds ±100kb clump regions from `04_gwas/1_mlma-loco`'s and `04_gwas/2_polmm`'s plink-clumped
GWAS output — one region set for dogCD (3 continuous factor traits from `1_mlma-loco`, plus all
14 categorical item traits from `2_polmm`, pooled together) and one for SIZE. Feeds
[`../cre_overlap/`](../cre_overlap) (Fig. 4e) — via two deposited reference copies
(`data/06_cCRE/dogCD_gws100kb_collapsed.bed` / `SIZE_gwas_clump100kb_regions_for_cCREintersect_clean.txt`),
not this process's own live output directly — see `cre_overlap/README.md`'s Notes for that caveat.

**2026-09-23**: extended to also pool in `2_polmm`'s 14 item traits (previously only the 3
mlma-loco factor traits + SIZE were included) — see Notes below.

**Pipeline engine:** Nextflow
**Containerized:** Yes (Seqera Wave) — `clump_regions` (r-base + bedtools)
**Estimated runtime:** [FILL IN]

Split out of the parent `cCRE.nf` (2026-09-15) into its own independent substage specifically so
`emission-state-plot`'s own stub test (deposited-inputs-only) wouldn't break the moment this
substage's `04_gwas`-dependent params became required — see the parent README and
`REPRODUCIBILITY_AUDIT.md`'s 2026-09-15 update.

---

## Input Data

No samplesheet — inputs are direct file params (see
[`nextflow_schema.json`](nextflow_schema.json)).

| Parameter | File | Source |
|-----------|------|--------|
| `clumping_dir` | `CCD{F1,F2,F3}_clump250.clumped` | `04_gwas/1_mlma-loco`'s published `clumping/` |
| `clumping_size_dir` | `SIZE_clump250.clumped` | `04_gwas/1_mlma-loco`'s published `clumping/size/` |
| `mlma_dir` | `CCD{F1,F2,F3}_..._LOCO.loco.mlma`, `SIZE_..._LOCO.loco.mlma` | `04_gwas/1_mlma-loco`'s published `gwas/` |
| `polmm_clumping_dir` | `polmmCCDitem<id>_clump250.clumped` (14 items) | `04_gwas/2_polmm`'s published `clumping/` |
| `polmm_gwas_dir` | `modi_simuMarkerOutput_POLMM_item<id>_FULLGRM_fullGeno.txt` (14 items) | `04_gwas/2_polmm`'s published `gwas/` |

---

## Output Data

| Path | Description |
|------|-------------|
| `regions/*` | ±100kb clump-region BED files, one per label (`dogCD`, `SIZE`) |
| `pipeline_info/versions.yml` | Per-process tool versions |

---

## Process Map

| # | Process | Description |
|---|---------|-------------|
| 1 | `BUILD_CLUMP_REGIONS_100KB` | Builds ±100kb regions around each clumped locus, per label (dogCD = F1+F2+F3+14 POLMM items pooled, SIZE = its own phenotype) |

---

## Notes

**POLMM's raw marker output has no `SNP` column.** `build_clump_regions_100kb.R` needs a
`SNP`+`bp` column pair to resolve secondary-clumped-SNP positions (plink's `.clumped` file only
lists them by ID, via `SP2`). `1_mlma-loco`'s `.loco.mlma` files have both natively; POLMM's
`modi_simuMarkerOutput_POLMM_item<id>_FULLGRM_fullGeno.txt` has `CHR`/`Position`/`bp` but no `SNP`
(GRAB.Marker never computes one). Fixed by constructing `SNP = CHR:Position` when missing —
matches the upstream `bcftools annotate --set-id '%CHROM:%POS'` convention that plink's own `SP2`
IDs are drawn from, confirmed this is what `CLUMP_ITEM` actually clumped against. Backward
compatible: the mlma-loco files already have a native `SNP` column, so this fallback never
triggers for them.

**Empirically validated against real data (2026-09-23).** No `04b-polmm` pixi-pipeline output
existed yet, but KatarinaTe pointed at real historical POLMM/mlma-loco output on Pelle
(`.../NEW_DA_IMP_MERGED/2024-10-11/`) — `build_clump_regions_100kb.R` was run directly (base R
only, no container needed) against the real `.clumped`/`modi_simuMarkerOutput...` files for all 3
factor traits and all 14 POLMM items, pooled together exactly as the live workflow now does.
Produced 84 real pooled "dogCD" regions, no errors — confirming the `SNP`-construction fallback
and the pooling logic both work correctly against real data, not just structurally. Locating the
real files required ruling out several decoy/duplicate candidates first (see
`REPRODUCIBILITY_AUDIT.md`'s 2026-09-23 update for the full trail — `_6PCs`/`.top1M.txt`
variants, and a `MODI_`-vs-`modi_` pair that turned out to be identical data, one with a
manually-added `SNP` column matching what the fallback constructs anyway).

**Genome-wide significance filter added (2026-09-23), and fully re-validated.** A `.clumped`
file only guarantees the lead SNP passed `--clump-p1` (1e-6, `CLUMP_FACTOR`/`CLUMP_ITEM`'s own
threshold — a loose "catch every candidate" cutoff, not genome-wide significance). Confirmed
directly with the manuscript author: the original GWAS-region reference files were built by
clumping broadly first, then manually selecting genome-wide significant clumps (lead SNP
P < 4e-7, the manuscript's own stated threshold) afterward. `build_clump_regions_100kb.R` now
automates that same filter — filter first, merge second. A "merge first, filter after"
alternative was also tested directly and produced a worse match, confirming filter-first is
correct (see `REPRODUCIBILITY_AUDIT.md`'s 2026-09-23 update).

Re-ran the real dogCD test above with the filter properly active: 84 → 36 real clumps (48 were
only suggestive, not genome-wide significant). Merged and compared against the deposited
`dogCD_gws100kb_collapsed.bed`: **0 bp missing** — every region in the deposited reference is
fully captured — plus exactly 400,000 bp of new coverage (two new 200kb regions,
`chr11:50174733-50374733` and `chr2:31264801-31464801`, contributed by POLMM items the old,
pre-POLMM deposited file never covered). Confirms the corrected pipeline reproduces the deposited
reference exactly, plus exactly the new coverage the POLMM-pooling scope change should add.

---

## Container

`clump_regions` — `r-base` + `bedtools` (the R script only uses base-R `read.table`/`write.table`;
`bedtools` does the final sort+merge). See
[`env/06_cCRE/README.md`](../../../env/06_cCRE/README.md).
