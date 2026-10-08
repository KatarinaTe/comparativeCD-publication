# 02. Data processing — 1. Survey data

Exploratory factor analysis and per-factor IRT fitting of the 14-item Darwin's Ark dogCD survey, producing per-dog F1/F2/F3 factor scores that [`3_Merging_filtering`](../3_Merging_filtering) merges with the genotype data.

**Pipeline engine:** Nextflow **Containerized:** Yes (Seqera Wave) **Estimated runtime:** \[FILL IN\]

------------------------------------------------------------------------

## Pipeline Overview

``` mermaid
%%{init: {'theme': 'neutral'}}%%
flowchart TD
    RESP(["response_df_dogCDitems.txt\n(14-item recoded checkpoint)"])
    QUEST(["Survey item metadata\n(questions.csv)"])

    EFA["RUN_EFA\n(diagnostics + fit,\nall reproduced, unused downstream)"]
    IRT["RUN_IRT_FACTOR\n(F1 / F2 / F3, one config-driven process)"]
    EXTRA["PLOT_EXTRA_ITEM_SE\n(F1 only)"]
    DIST["PLOT_FACTOR_DISTRIBUTIONS"]

    RESP --> EFA
    QUEST --> EFA
    RESP --> IRT
    IRT --> DIST
    IRT -->|F1 only| EXTRA

    OUT[("F1/F2/F3_CCD3F.txt")]
    IRT --> OUT
    NEXT(["3_Merging_filtering"])
    OUT --> NEXT
```

------------------------------------------------------------------------

## Input Data

No samplesheet — inputs are direct file params (see [`nextflow_schema.json`](nextflow_schema.json)).

| Parameter | File | Source |
|-------------------------|------------------|-----------------------|
| `response_df` | `response_df_dogCDitems.txt` | Download from [doi.org/10.17044/scilifelab.33339309](https://doi.org/10.17044/scilifelab.33339309) into `data/02_data_processing/` — already-recoded 14-item checkpoint (25,302 dogs) |
| `survey_questions_csv` | `DarwinsArk_20220315_questions.csv` | [`data/02_data_processing/DarwinsArk_20220315_questions.csv`](../../../data/02_data_processing/DarwinsArk_20220315_questions.csv) |

------------------------------------------------------------------------

## Output Data

| Path | Description |
|------------------------------|-------------------------------------------------|
| `scores/F{1,2,3}_CCD3F.txt` | Per-dog factor score + SE + item responses, one file per factor |
| `efa/diagnostics/efa_diagnostics.log` | nfactors/KMO/Bartlett/determinant/cronbach's alpha/loadings — reproduced, feeds nothing downstream |
| `efa/plots/*.pdf` | Parallel-analysis, nScree, and uniqueness plots |
| `irt/diagnostics/F{1,2,3}_diagnostics.log` | Per-factor mirt/CTT diagnostics |
| `irt/plots/*.pdf` | Per-factor trace/item-info/test-info/EAP-MAP/NA-count plots (+ F1's extra plot) |
| `distributions/diagnostics/distribution_diagnostics.log` | Per-factor score/SD ratio diagnostics |
| `distributions/plots/factor_distribution_F{1,2,3}.pdf` | Final per-factor score histograms |
| `pipeline_info/versions.yml` | Per-process tool versions |

------------------------------------------------------------------------

## Process Map

| \# | Process | Description |
|----------------|-----------------------|----------------------------------|
| 1 | `RUN_EFA` | nfactors/parallel-analysis/nScree diagnostics, KMO/Bartlett/determinant/cronbach's alpha, EFA fit + uniqueness plot |
| 2 | `RUN_IRT_FACTOR` | Per-factor (F1/F2/F3) IRT fit, CTT solution, factor scores, and the factor's saved CSV |
| 3 | `PLOT_EXTRA_ITEM_SE` | F1-only diagnostic scatter plot (an item's response vs. factor SE); split into its own process because Nextflow output arity can't express "produced only for some factors" as an optional output on a process all three factors share. Reads F1's already-saved CCD3F.txt rather than recomputing |
| 4 | `PLOT_FACTOR_DISTRIBUTIONS` | Final per-factor score histograms (barrier — needs all three factors) |

------------------------------------------------------------------------

## Notes

- `RUN_EFA`'s diagnostics and each factor's IRT diagnostics feed nothing downstream — they are captured as real log/plot files rather than left to an interactive session.
- `RUN_IRT_FACTOR` is one config-driven process handling all three factors, since the per-factor differences are genuine and are passed through as config: `M2()`'s `type` (`"M2*"` for F1/F3, `"C2"` for F2), which `fscores()` method feeds the saved CSV (MAP for F1/F3, EAP for F2 — both EAP and MAP are always computed and plotted for every factor, only the one joined into the final output differs), and per-factor `theta_range`/color-palette values for the trace/test-info/item-info plots. F1's one extra diagnostic scatter plot (`item7` vs. factor SE) is its own process, `PLOT_EXTRA_ITEM_SE`, invoked only for F1 — see the Process Map above for why.
- Per-factor items are selected by column from the already-recoded `response_df_dogCDitems.txt` checkpoint rather than re-derived from raw survey data: each factor's recoding only touches that factor's own items, so this produces the same values as re-deriving from raw data.
- `nFactors::parallel()`'s random simulation (`rep=100`, feeding the nScree diagnostic plot) has no set seed and is not reproducible run-to-run. This affects only that one diagnostic plot — the number of factors used (`nFac = 3`) is a fixed value, and the EFA fit itself is a deterministic iterative fit given fixed inputs.
- 3 factors are used, not 4: the 4th factor's SS loadings (0.669) fall below 1, the standard retention threshold; the three used factors all have SS loadings above 1 (2.126, 1.357, 1.063).
- Each factor's `item.weights` vector is a hardcoded literal (`FACTOR_CONFIGS` in the main workflow), not recomputed from `mirt()`'s fitted parameters at runtime.
- `M2()` errors for every factor (too few degrees of freedom); the error is caught and logged, and execution continues to `itemfit()` and beyond.
- `run_irt_factor.R` needs an explicit `library(ggplot2)` call, since this pipeline only loads the packages actually used past the `response_df_dogCDitems.txt` checkpoint (see `env/02_data_processing/README.md`).
- A handful of small per-factor asymmetries in diagnostic-log print ordering and in the EAP/MAP plot's exact per-factor title/xlim/color are not reproduced (none affect the computed scores or the saved `CCD3F.txt`); the EAP/MAP and NA-count plots use one consistent generic form across all three factors.

------------------------------------------------------------------------

## Container

Only the R packages actually used past the `response_df_dogCDitems.txt` checkpoint are included — see [`env/02_data_processing/README.md`](../../../env/02_data_processing/README.md) for the full package list and how `nFactors`/`ltm`/`ggmirt` (none of which are on conda-forge) are layered on top of the conda-resolved base via Wave's `--conda-run-command`.