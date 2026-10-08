# Container Derivation for 04 GWAS Processes

This directory contains Conda environment specifications used to freeze Docker containers for the
Nextflow processes defined under [`code/04_gwas/`](../../code/04_gwas/).

## Container Types

As in [`env/01_mapping`](../01_mapping/README.md), containers here are frozen with the
[Seqera Wave CLI](https://github.com/seqeralabs/wave-cli), which turns a Conda environment spec
into a reproducible Docker image.

### Commands

```bash
~/Downloads/wave-1.8.2-macos-arm64 --conda-file env/04_gwas/2_polmm/polmm_grab.yml --platform linux/amd64 --freeze --await
```

Frozen to `community.wave.seqera.io/library/polmm_grab:4c906c4feeeb8621`, substituted into
`RUN_POLMM_ITEM` in `code/04_gwas/2_polmm/process_definitions.nf`.

```bash
/tmp/wave-cli/wave-1.8.2-macos-arm64 --conda-file env/04_gwas/3_finemap_susie/polyfun.yml \
  --conda-run-command "RUN git clone https://github.com/omerwe/polyfun.git /opt/polyfun && cd /opt/polyfun && git checkout 3e657a1066" \
  --conda-run-command "RUN chmod +x /opt/polyfun/*.py" \
  --platform linux/amd64 --freeze --await -o json
```

Frozen to `community.wave.seqera.io/library/polyfun:cb5ec2572b0e5a72`, substituted into
`HARMONIZE_SUMSTATS`/`MUNGE_SUMSTATS`/`FINEMAP_REGION` in
`code/04_gwas/3_finemap_susie/process_definitions.nf`.

`04_gwas/1_mlma-loco` needs no new container — it reuses `3_Merging_filtering`'s
`merging_filtering_tools` (has `gcta==1.94.1`) for mlma-loco/sumstats-export, and the existing
`plink` container (derivation: [`env/shared/README.md`](../shared/README.md)) for clumping.
`2_polmm`'s own clumping/sumstats-export steps do the same, and `3_finemap_susie` reuses
`merging_filtering_tools` (has R) for its genofile-creation R script and the same shared `plink`
container for its own plink filtering steps.

## Environment Files

### [`2_polmm/polmm_grab.yml`](./2_polmm/polmm_grab.yml)

Single-purpose container for `RUN_POLMM_ITEM` — the `GRAB` R package's `GRAB.NullModel`/
`GRAB.Marker` functions (POLMM method, for ordinal per-item GWAS).

| Tool | Version | Purpose |
|------|---------|---------|
| r-base | 4.5.3 | Runtime |
| r-grab | 0.2.5 | `GRAB.NullModel`/`GRAB.Marker` (POLMM ordinal-trait GWAS) |

Validated with `pixi lock --dry-run` against `conda-forge`/`bioconda` for `linux-64` — resolves
cleanly. Confirmed both `GRAB.NullModel` and `GRAB.Marker` load and exist under those exact names
in the built container.

GRAB is pinned to 0.2.5, the newest version available on conda-forge (back to 0.2.2) and CRAN's
archive (back to 0.2.1). The version used in the original analysis may have been earlier (possibly
0.1.1); its GitHub repository could not be located to check algorithm details against that
possibility.

### [`3_finemap_susie/polyfun.yml`](./3_finemap_susie/polyfun.yml)

Container for `HARMONIZE_SUMSTATS`/`MUNGE_SUMSTATS`/`FINEMAP_REGION` — [PolyFun](https://github.com/omerwe/polyfun)
(`munge_polyfun_sumstats.py`, `finemapper.py`, SuSiE via `r-susier`/`rpy2`). PolyFun has no
versioned releases (rolling git history only) and isn't pip/conda-installable as a package, so its
own `polyfun.yml` conda spec is layered with a `git clone` + `git checkout` build step.

| Tool | Version | Purpose |
|------|---------|---------|
| python | 3.8 | Runtime (PolyFun's own pin) |
| r-susier | 0.11.92 | SuSiE finemapping, called from `finemapper.py` via `rpy2` |
| pandas / pyarrow / scipy / scikit-learn / pandas-plink / bgen | PolyFun's own pins | sumstats munging, genotype reading |
| openpyxl | latest on conda-forge | reads `data/04_gwas/SIZE_260603_CSnonfunct.xlsx` in `extract_size_regions.py` (added 2026-09-02, not part of PolyFun's own env) |
| PolyFun (git) | commit `3e657a1066` | `munge_polyfun_sumstats.py`, `finemapper.py` (cloned to `/opt/polyfun`) |

Rebuilt 2026-09-02 (tag `cb5ec2572b0e5a72`, was `4ed13d7a6dd05749`) solely to add `openpyxl` —
also substituted into `4_qc_signals_IGV`'s processes, which reuse this same container.

Validated with `pixi lock --dry-run` against `conda-forge` for `linux-64` — resolves cleanly.
Confirmed in the built container: PolyFun's scripts are present at `/opt/polyfun`, at the pinned
commit; `rpy2` successfully imports and loads the `susieR` R package (the exact bridge
`finemapper.py --method susie` depends on); `openpyxl` imports successfully.

PolyFun has no tags/releases, so it is pinned to commit `3e657a1066` (2024-07-28) — the latest
commit that existed before the original analysis's own recorded run date (its working directory is
named `.../2024-10-11/`). This is a traceable anchor but not a confirmed exact match; an earlier
commit may have been used.

## Notes

- All containers are built for `linux-64` platform.
- Channels: `conda-forge` and `bioconda`.
