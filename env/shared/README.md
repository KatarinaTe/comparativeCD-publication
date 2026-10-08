# Container Derivation for Containers Shared Across Stages

Some single-tool containers are referenced by more than one stage's `process_definitions.nf`
rather than living under one stage's own `env/<stage>/`. This directory is their one home, so the
build derivation is documented once instead of repeated per consumer.

## [`plink.yml`](./plink.yml)

`bioconda::plink=1.90b7.7`. Used by `01_mapping/2_cohort_genotyping`, `02_data_processing/
3_Merging_filtering`, and every `04_gwas` substage.

## [`glimpse-bio.yml`](./glimpse-bio.yml)

`bioconda::glimpse-bio=1.1.1`. Used by `01_mapping/2_cohort_genotyping` and `03_qc`.

## Provenance and how to rebuild (`plink`/`glimpse-bio`)

Both were originally built via the Seqera Containers web interface, not the Wave CLI directly, so
no build command was ever recorded at the time. Rather than leave that gap, the exact spec and
resolved lock file were recovered straight from the running images themselves — every container
Wave builds this way bakes in its own input spec and resolved package list:

```bash
cid=$(docker create --platform linux/amd64 community.wave.seqera.io/library/plink:1.90b7.7--feb905bcf3811778)
docker cp "$cid:/tmp/conda.yml" env/shared/plink.yml          # exact input spec Wave used
docker cp "$cid:/tmp/environment.lock" env/shared/plink.lock  # exact resolved package list
docker rm "$cid"
```

(same for `glimpse-bio`, substituting its image reference,
`community.wave.seqera.io/library/glimpse-bio:1.1.1--7517c3cd783f15ab`). `plink.yml`/`glimpse-bio.yml` above are
exactly what came out of that — not reconstructed from guesswork.

Confirmed via `wave --inspect` against both tags: both share the same `mambaorg/micromamba` base
layers as every other container in this project, i.e. genuine single-package Wave builds, not an
nf-core/biocontainers passthrough under a Wave-shaped name.

**Verified, and important: re-running `wave --conda-file env/shared/plink.yml --platform linux/amd64
--freeze --await` today does *not* reproduce the pinned tag above, or even the same environment.**
It resolves successfully, but to a different image (a fresh build was actually run to check this)
with newer transitive dependencies — `libgcc-16.1.0`/`libopenblas-0.3.34` instead of the original
build's `libgcc-14.2.0`/`libopenblas-0.3.29`, etc. `plink=1.90b7.7` itself is the one thing pinned
exactly in the `.yml`; everything underneath it (compiler runtime, BLAS, etc.) floats to "latest
compatible with conda-forge's index at build time" unless also pinned — and conda-forge's index has
moved on since the original build. **The `.lock` file, not the `.yml`, is what actually guarantees
the original exact environment**: it's the literal resolved package list (exact build strings,
exact URLs) from the real historical image, recoverable offline and forever, independent of what
conda-forge's index looks like whenever someone reads this:

```bash
micromamba create -n plink -f env/shared/plink.lock   # or: conda create --file env/shared/plink.lock
```

This installs the exact package set the original container had — not a fresh Wave build, but a
locally-reconstructable environment with the identical resolved versions, which is the strongest
reproducibility guarantee available here.

**One build-process quirk worth knowing**: both `.lock` files include a `procps-ng` package that
isn't in the input `.yml` at all — confirmed via each image's own `/opt/conda/conda-meta/history`,
this is Wave's own build template installing it as a second step after resolving the requested
spec, not something a conda file needs to request. Expect it in the resolved lock; don't try to
add it to the input spec.

## [`gwas_supplement_plots.yml`](./gwas_supplement_plots.yml)

Container for `PLOT_GWAS_RESULTS` (`04_gwas/5_plotting`), `PLOT_HERITABILITY`
(`07_additional_plots_and_analyses/1_heritability_plot`),
`REPRODUCE_GWAS_CATALOG_OVERLAP` (`07_additional_plots_and_analyses/2_gwas_catalog_overlap`), and
`PLOT_NETCOLOC_SCATTER` (`05_network/1_NetColoc`, Fig. 1g — added 2026-09-21):
`r-base`, `r-tidyverse`, `r-data.table`, `r-readxl`, `r-ggrepel`, `r-qqman` — all conda-forge/
bioconda, no exotic dependencies. Reused across all 4 rather than building separate containers,
same reasoning `02_data_processing/4_Plotting` used to reuse `3_Merging_filtering`'s container.
(`1_NetColoc`'s own `netcoloc` container is pure Python — no R at all — so this R-only plotting
step reuses this one instead, rather than adding R to `netcoloc` or building a 5th container.)

Precheck (pixi dry-run solver, no install/build) confirmed the spec resolves cleanly before the
real build.

Built with the [Seqera Wave CLI](https://github.com/seqeralabs/wave-cli) (v1.8.2):

```bash
wave --conda-file env/shared/gwas_supplement_plots.yml --platform linux/amd64 --freeze --await -o json
```

Frozen to `community.wave.seqera.io/library/gwas_supplement_plots:e626727dcb23673b`.

**No `.lock` file for this one** — checked directly (`docker run ... find / -iname
'*environment*'`), confirmed empty. Unlike `plink`/`glimpse-bio` above, this is a **multi-package**
Wave `--freeze` build (6 top-level packages, ~100 resolved), and that path does not bake a
`/tmp/conda.yml`/`/tmp/environment.lock` into the image the way a single-package build does —
verified against `go_semantic_plots` too (same empty result). `AGENTS.md`'s convention describing
these files as baked into "every Wave-built container" should be narrowed to single-package
builds specifically; multi-package `--freeze` builds like this one and `go_semantic_plots` have
no such file to recover, so there's no lock file to add here, not an oversight.

## Notes

- Built for `linux-64` platform.
- Channels: `conda-forge`, `bioconda`.
