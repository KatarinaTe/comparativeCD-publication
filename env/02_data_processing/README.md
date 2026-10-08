# Container Derivation for 02 Data Processing Processes

This directory contains Conda environment specifications used to freeze Docker containers for the
Nextflow processes defined under [`code/02_data_processing/`](../../code/02_data_processing/).

## Container Types

As in [`env/01_mapping`](../01_mapping/README.md), containers here are frozen with the
[Seqera Wave CLI](https://github.com/seqeralabs/wave-cli), which turns a Conda environment spec
into a reproducible Docker image.

### Commands

```bash
~/Downloads/wave-1.8.2-macos-arm64 --conda-file env/02_data_processing/2_Axiom_imputation/axiom_liftover_impute.yml --freeze --await
```

Frozen to `community.wave.seqera.io/library/axiom_liftover_impute:27db68ae33d39b0e`, already
substituted into every process in `code/02_data_processing/2_Axiom_imputation/process_definitions.nf`.
Re-run the command above (bump the Wave CLI version if a newer release exists) and update the tag
in `process_definitions.nf` if `axiom_liftover_impute.yml` ever changes.

```bash
~/Downloads/wave-1.8.2-macos-arm64 \
  --conda-file env/02_data_processing/1_Survey_data/survey_data_efa_irt.yml \
  --conda-run-command "RUN Rscript -e \"install.packages(c('nFactors','ltm'), repos='https://cloud.r-project.org')\"" \
  --conda-run-command "RUN Rscript -e \"remotes::install_github('masurp/ggmirt')\"" \
  --platform linux/amd64 --freeze --await
```

Frozen to `community.wave.seqera.io/library/survey_data_efa_irt:f4e0f4b9755aa129`, already
substituted into every process in `code/02_data_processing/1_Survey_data/process_definitions.nf`.
The `--conda-run-command` steps layer on the two packages that have no conda-forge build at all
(see below) — `--platform linux/amd64` is passed explicitly since this one was built from a
macOS ARM64 dev machine.

```bash
~/Downloads/wave-1.8.2-macos-arm64 --conda-file env/02_data_processing/3_Merging_filtering/merging_filtering_tools.yml --platform linux/amd64 --freeze --await
```

Frozen to `community.wave.seqera.io/library/merging_filtering_tools:0e56b6d8b7fd7863`, substituted
into every process in `code/02_data_processing/3_Merging_filtering/process_definitions.nf` that
needs R/KING/GCTA (plink-only processes use a separate, existing `plink` container instead —
derivation: [`env/shared/README.md`](../shared/README.md)). Also
reused as-is by `4_Plotting` — same R/dplyr/ggplot2/patchwork/ggtext stack, no separate build
needed.

## Environment Files

### [`2_Axiom_imputation/axiom_liftover_impute.yml`](./2_Axiom_imputation/axiom_liftover_impute.yml)

Multi-tool container for the whole `2_Axiom_imputation` pipeline. Includes:

| Tool | Version | Purpose |
|------|---------|---------|
| plink | 1.90b7.7 | canFam3 filtering, liftover marker bookkeeping, final PLINK conversion/merge |
| plink2 | 2.0.0a.6.9 | Per-chromosome VCF recoding with reference-based allele orientation |
| ucsc-liftover | 482 | canFam3 -> canFam4 coordinate liftover |
| bcftools | 1.20 | VCF/BCF manipulation |
| htslib | 1.20 | HTS file I/O |
| openjdk | 17.0.18 | Runs the externally-fetched `conform-gt`/`beagle` jars (see below) |
| gawk | 5.4.1 | Marker bookkeeping file derivation |

Validated with `pixi lock --dry-run` against `conda-forge`/`bioconda` for `linux-64` — resolves
cleanly with all versions pinned exactly as listed above.

`conform-gt` and `beagle` are deliberately **not** included in this environment. `conform-gt` has
no bioconda package at all, and while bioconda does package `beagle`, it currently resolves to a
2025 build (`>=5.5_27feb25.75f`) rather than the `22Jul22.46e` release the original pipeline
actually ran — a different beagle version can change imputation results, not just repackage the
same one. Both jars are fetched instead by `pixi run 02a-fetch-axiom-jars` into
`code/02_data_processing/2_Axiom_imputation/bin/` (gitignored) and passed into the relevant
processes as regular file inputs.

### [`1_Survey_data/survey_data_efa_irt.yml`](./1_Survey_data/survey_data_efa_irt.yml)

Container for the whole `1_Survey_data` EFA/IRT pipeline. Only channel-forge (`conda-forge`)
packages are declared in the yml itself — `nFactors`, `ltm`, and `ggmirt` have no conda-forge
build (and `nFactors`/`ltm` exist only on Anaconda's `r` channel, which falls under the same
commercial-license restriction as `defaults` — not used here) — so they're installed from CRAN/
GitHub via `--conda-run-command` on top of the conda-resolved base, layered in the order:
`r-remotes` (from conda) -> `install.packages(c('nFactors','ltm'))` -> `remotes::install_github`.

| Tool | Version | Source | Purpose |
|------|---------|--------|---------|
| r-base | 4.5.3 | conda-forge | Runtime |
| dplyr | 1.2.1 | conda-forge | `%>%`, `mutate`/`filter`/`left_join`/`arrange`/`select` |
| data.table | 1.18.4 | conda-forge | `data.table()`, `setnames()` |
| psych | 2.6.5 | conda-forge | `KMO`, `cortest.bartlett`, `fa` |
| mirt | 1.47 | conda-forge | IRT model fitting, `fscores` |
| lavaan | 0.7-2 | conda-forge | CTT `cfa`/`standardizedsolution` |
| ggplot2 | 4.0.3 | conda-forge | Distribution plots |
| remotes | 2.5.0 | conda-forge | Installs `ggmirt` from GitHub |
| nFactors | (CRAN latest) | CRAN, via `install.packages()` | Parallel analysis, `nScree` |
| ltm | (CRAN latest) | CRAN, via `install.packages()` | `cronbach.alpha` |
| ggmirt | (GitHub `masurp/ggmirt`, latest) | `remotes::install_github()` | `tracePlot`/`itemInfoPlot`/`testInfoPlot` |

Only the packages actually used past the `response_df_dogCDitems.txt` checkpoint are included —
the original script's own library header (`phateR`, `ggpubr`, `zoo`, `stringr`, `tm`, `caret`,
`FactoMineR`, `mice`, etc.) loads a much larger stack, but almost all of it is only used in the
raw-survey-CSV ingestion this pipeline doesn't reproduce (see the pipeline README).

### [`3_Merging_filtering/merging_filtering_tools.yml`](./3_Merging_filtering/merging_filtering_tools.yml)

Container for the R/KING/GCTA portions of `3_Merging_filtering` (and all of `4_Plotting`, which
reuses it directly).

| Tool | Version | Purpose |
|------|---------|---------|
| king | 2.3.2 | Relatedness/duplicate detection |
| gcta | 1.94.1 | Genetic relationship matrix (GRM) |
| r-base | 4.5.3 | Runtime |
| r-dplyr | 1.2.1 | `left_join`/`full_join`/`mutate`/`filter` |
| r-ggplot2 | 4.0.3 | Distribution/histogram/PCA plots |
| r-patchwork | 1.3.2 | Combining multiple ggplot panels into one figure |
| r-ggtext | 0.1.2 | `element_textbox_simple()` styled plot titles/axis labels |

## Notes

- All containers are built for `linux-64` platform.
- Channels: `conda-forge` and `bioconda`.
