# Container Derivation for 06_cCRE Processes

This directory contains the Conda environment specification used to freeze the Docker container for
the Nextflow process defined under [`code/06_cCRE/`](../../code/06_cCRE/).

## [`emission_state.yml`](./emission_state.yml)

Container for `EMISSION_STATE_PLOT`, plus `cre_overlap/`'s `PLOT_TISSUE_SPECIFICITY` and
`PLOT_REGION_FOREST` (added 2026-09-21 — both are simple ggplot2/dplyr forest plots, no new
packages needed): `r-base`, `r-ggplot2`, `r-patchwork`, `r-dplyr`, `r-tidyr`, `r-scales` — all
conda-forge, no exotic dependencies.

Built with the [Seqera Wave CLI](https://github.com/seqeralabs/wave-cli):

```bash
wave --conda-file env/06_cCRE/emission_state.yml --platform linux/amd64 --freeze --await -o json
```

Frozen to `community.wave.seqera.io/library/emission_state:0737d6956e5a4a7d`.

### Why R, not the duplicate Python notebook

`emission_state_UU.R` and `emission_state_UU.ipynb` are two independent implementations of the same
figure, from the same 2 deposited inputs. The R version was chosen after confirming it's what
actually produced the deposited historical output: the deposited `emission_states_plot.png` is
4800x2100px, exactly matching the R script's own `ggsave(width=16, height=7, dpi=300)` call — the
Python `savefig` call produces different pixel dimensions entirely. Validated by running inside the
built container and comparing pixel-for-pixel against the deposited PNG: identical dimensions,
colors, layout, and all 9 chromatin states — the only difference is minor tick-label spacing on the
"Genome %" colorbar, a font-availability rendering quirk between R installations, not a logic bug.

## [`clump_regions.yml`](./clump_regions.yml)

Container for `BUILD_CLUMP_REGIONS_100KB`, plus `cre_overlap/`'s `BUILD_EPICDOG_TISSUE_OVERLAPS`
and `BUILD_UU_REGION_BP_OVERLAPS` (added 2026-09-21 — both are pure bedtools/awk shell scripts):
`r-base` + `bedtools` — no other dependencies (the R script only uses
`read.table`/`write.table`, no packages beyond base R; `bedtools` does the final sort+merge).

Built with the [Seqera Wave CLI](https://github.com/seqeralabs/wave-cli) (v1.8.2):

```bash
wave --conda-file env/06_cCRE/clump_regions.yml --platform linux/amd64 --freeze --await -o json
```

Frozen to `community.wave.seqera.io/library/clump_regions:09876be1daa6bcb4`.

## [`compare_UU_epicdog.yml`](./compare_UU_epicdog.yml)

Container for `COMPARE_UU_EPICDOG_ELEMENTS` (Fig. 4d, added 2026-09-22): `r-base`, `r-readr`,
`r-dplyr`, `r-tidyr`, `r-ggplot2`, `r-valr` — all conda-forge (`r-valr` doesn't need `bioconda`
here). The source script's `library(tidyverse)` was swapped for these four specific packages in
the wired `bin/` copy — no need for the full `tidyverse` meta-package.

Built with the [Seqera Wave CLI](https://github.com/seqeralabs/wave-cli):

```bash
wave --conda-file env/06_cCRE/compare_UU_epicdog.yml --platform linux/amd64 --freeze --await -o json
```

Frozen to `community.wave.seqera.io/library/compare_uu_epicdog:ebff580a4fa19c19`.

## Notes

- Built for `linux-64` platform.
- Channels: `conda-forge`, `bioconda` (the latter only needed for `clump_regions.yml`'s
  `bedtools`).
