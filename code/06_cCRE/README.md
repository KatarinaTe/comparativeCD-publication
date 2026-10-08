# 06. cCRE

Candidate cis-regulatory element (cCRE) annotation and dogCD/EpicDog comparisons behind Fig. 4c–f.
A small, real Nextflow pipeline (4 independent substages) sits alongside a handful of loose,
paused, or superseded scripts. This README is deliberately honest about which is which rather
than presenting it as more finished than it is.

**Pipeline engine:** Nextflow (for the 4 substages below)
**Containerized:** Yes, for the Nextflow substages
**Estimated runtime:** [FILL IN]

---

## Live Nextflow pipeline (4 independent substages)

| Substage | Process(es) | Fig. | Description |
|---|---|---|---|
| This directory (`cCRE.nf`) | `EMISSION_STATE_PLOT` | 4c | 9-state emission-probability heatmap from deposited-only inputs |
| [`clump_regions_100kb/`](clump_regions_100kb) | `BUILD_CLUMP_REGIONS_100KB` | — (feeds `cre_overlap/` below) | Builds ±100kb clump regions from `04_gwas/1_mlma-loco`'s and `04_gwas/2_polmm`'s output |
| [`cre_overlap/`](cre_overlap) | 4 processes (2 bedtools bp-overlap builders, 2 R plotters) | 4e, 4f | GWAS-region overlap with cCREs: EpicDog tissue-specificity (4e) and UU per-region (4f) |
| [`compare_UU_epicdog/`](compare_UU_epicdog) | `BUILD_ALL_TISSUES_FILTERED`, `COMPARE_UU_EPICDOG_ELEMENTS` | 4d | UU vs. EpicDog element-count comparison (Shared/Unique, by functional class) |

Split into independent substages specifically because each has a different input-availability
profile: `emission-state-plot`'s inputs are deposited-only (no pipeline dependency),
`clump-regions-100kb` needs `1_mlma-loco`'s real output, `cre_overlap` needs both fetched
chromatin-state BEDs and the GWAS-region files, and `compare_UU_epicdog` needs its own live
EpicDog fetch (`06d-fetch-epicdog-bed`) plus a deposited UU-side reference copy — bundling any of
these together broke another's own stub test the moment an unrelated substage's params became
mandatory. See `REPRODUCIBILITY_AUDIT.md`'s 2026-09-15 and 2026-09-21 updates for the full story.

```bash
pixi run 06a-emission-state-plot
pixi run 06b-clump-regions-100kb   # needs 04a-mlma-loco + 04b-polmm
pixi run 06c-cre-overlap           # needs 06b-clump-regions-100kb
pixi run 06f-compare-uu-epicdog    # needs 06d-fetch-epicdog-bed + 06e-fetch-chromatin-states-bed
                                    # — both pulled in automatically
```

---

## Loose scripts — status of each

Fig. 4d's `compare_UU_epicdog_elements.R` was folded into the Nextflow pattern above on
2026-09-22 (see `compare_UU_epicdog/`) — the original stays here, untouched, for provenance; the
wired copy lives at `compare_UU_epicdog/bin/`. What each remaining loose script actually is:

| Script | Status | Role |
|---|---|---|
| `compare_UU_epicdog_elements.R` | **Wired into Nextflow (2026-09-22)**, original kept here untouched | Fig. 4d comparison script (UU vs. EpicDog element counts by Promoters/Enhancers/Repressed). Two real bugs found getting it to actually run — a `valr::read_bed()` parsing crash and its `n_fields` defaulting to 3 (so `filter(name %in% ...)` errored, column 4 came back as `X4` not `name`) — see `compare_UU_epicdog/bin/`'s header and README. Output confirmed byte-for-byte identical to the already-deposited `data/06_cCRE/EPIC_comparison.txt`. |
| `build_all_tissues_filtered.sh` | **Wired into Nextflow (2026-09-23)**, original kept here untouched | Reconstructs `data/06_cCRE/all_tissues_filtered.bed` (the UU-side input `compare_UU_epicdog/` now builds live via `BUILD_ALL_TISSUES_FILTERED`, no longer a static deposited copy) from the 8 fetched per-region UU BEDs — live-built version has 7 fewer lines than the old deposited copy (stray `track name=...` header artifacts, confirmed analytically inert either way). |
| `compare_UU_epicdog_elements_bedtools.sh`, `diagnose_other_tissues.sh` | **Archived (2026-09-22)** | One-off cross-checks written while finding the mammary-gland bug and validating Fig. 4d's logic against an independent tool. Never part of any pipeline. Moved to [`archive/cCRE_fig4_debugging_aids/`](../../archive/cCRE_fig4_debugging_aids) — see `REPRODUCIBILITY_AUDIT.md`'s 2026-09-17 and 2026-09-22 updates. |
| `testing_CRE_overlap.R`, `testing_CRE_overlap_260617.R`, `check_CRE_overlap.sh`, `SIZE_check_UU_CRE_overlap.sh` | **Superseded — archived (2026-09-21)** | The real cCRE-overlap logic these contained (Fig. 4e and the EpicDog tissue-specificity odds ratios) is now in [`cre_overlap/`](cre_overlap) — see its README's Notes for exactly what was kept vs. dropped (most of `testing_CRE_overlap_260617.R` was exploratory/superseded duplicate work, never part of either published figure). Moved to [`archive/testing_CRE_overlap_superseded/`](../../archive/testing_CRE_overlap_superseded). |
| `emission_state_UU.R` / `.ipynb` | **`emission_state_UU.R` is the real source, `.ipynb` is a superseded duplicate** | Confirmed by matching the deposited `emission_states_plot.png`'s exact pixel dimensions against each script's own `ggsave()`/`savefig()` call — see `env/06_cCRE/README.md`. `.R` is what `EMISSION_STATE_PLOT` actually runs. |

---

## Known gaps

- `cre_overlap/` reads the same 8 per-region UU BEDs `compare_UU_epicdog/` fetches live via
  `06e-fetch-chromatin-states-bed`, but isn't (yet) wired to that same task — its `uu_states_dir`
  param still points at the deposited copy directly, a follow-up not done here.

---

## Container

See [`env/06_cCRE/README.md`](../../env/06_cCRE/README.md) for the frozen containers
(`emission_state`, `clump_regions`, `compare_uu_epicdog`) and their derivation — `cre_overlap/`
reuses the first two, no new build; `compare_UU_epicdog/` has its own.
