# Container Derivation for 07_additional_plots_and_analyses/3_cross_species_expression

## [`cross_species_expression.yml`](./cross_species_expression.yml)

Container for `CROSS_SPECIES_EXPRESSION`: `r-base`, `r-seurat`, `r-dplyr`, `bioconductor-fgsea`,
`bioconductor-ucell`, `bioconductor-complexheatmap`, `r-circlize`. The one genuinely fragile
container in the `07_additional_plots_and_analyses` conversion batch — Seurat's own dependency
tree resolves to ~150 transitive packages.

Precheck (pixi dry-run solver, no install/build) confirmed the full spec resolves cleanly,
including all 4 non-trivial packages (`bioconductor-fgsea`, `bioconductor-ucell`,
`bioconductor-complexheatmap`, `r-circlize`), before the real build.

Built with the [Seqera Wave CLI](https://github.com/seqeralabs/wave-cli) (v1.8.2):

```bash
wave --conda-file env/07_additional_plots_and_analyses/3_cross_species_expression/cross_species_expression.yml --platform linux/amd64 --freeze --await -o json
```

Frozen to `community.wave.seqera.io/library/cross_species_expression:3b296196c0c71539`.

No `.lock` file — same reason as `env/shared/gwas_supplement_plots.yml` (see that README): a
real multi-package Wave `--freeze` build doesn't bake `/tmp/conda.yml`/`/tmp/environment.lock`
into the image, checked directly.

**Pinned 2026-09-22** (previously unpinned — `.yml` had no version numbers at all, so a rebuild
from this spec today could silently resolve to different versions than what's actually frozen in
the already-built image). Recovered the exact resolved versions directly from the running
container (`Rscript -e 'packageVersion(...)'` for each package) rather than guessing:
`r-base==4.5.3`, `r-seurat==5.5.1`, `r-dplyr==1.2.1`, `bioconductor-fgsea==1.36.2`,
`bioconductor-ucell==2.14.0`, `bioconductor-complexheatmap==2.26.1`, `r-circlize==0.4.18`.
Re-verified the fully-pinned spec still resolves cleanly via `pixi lock --dry-run` before
committing it — same versions, confirming this doesn't require a rebuild, just makes the `.yml`
accurately describe what's already in the frozen image.

Note: the manuscript's own Methods text cites `Seurat v5.3.0`/`ComplexHeatmap v2.20.0` for this
analysis (and separately `Seurat v4` for a related one) — close to, but not exactly, what's
pinned here (`5.5.1`/`2.26.1`). Those specific historical versions weren't available when this
container was built; see `REPRODUCIBILITY_AUDIT.md`'s 2026-09-22 update for the broader
tool-version reconciliation this is part of.

## Notes

- Built for `linux-64` platform.
- Channels: `conda-forge`, `bioconda`.
