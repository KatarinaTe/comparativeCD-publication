# Container Derivation for 05 Network Processes

This directory contains Conda environment specifications used to freeze Docker containers for the
Nextflow processes defined under [`code/05_network/`](../../code/05_network/).

## Container Types

As in [`env/04_gwas`](../04_gwas/README.md), containers here are frozen with the
[Seqera Wave CLI](https://github.com/seqeralabs/wave-cli), which turns a Conda environment spec
into a reproducible Docker image.

### Commands

```bash
wave --conda-file env/05_network/1_NetColoc/netcoloc.yml \
  --conda-run-command "RUN pip install --no-deps ndex-dev" \
  --conda-run-command "RUN git clone https://github.com/michaelkyu/ddot.git /opt/ddot && cd /opt/ddot && git checkout v1.0" \
  --conda-run-command "RUN pip install --no-deps /opt/ddot" \
  --platform linux/amd64 --freeze --await -o json
```

Frozen to `community.wave.seqera.io/library/netcoloc:8ec66694a6aee19e`.

This one container serves **both** `1_NetColoc` and `2_CrossSpeciesBMI` — not scope creep: the
package pins below come directly from `CrossSpeciesBMI`'s own `environment.yml`
(https://github.com/sarah-n-wright/CrossSpeciesBMI), which is a full pip-freeze lockfile that
includes packages (`ddot`, `cdapsutil`, `gprofiler-official`) only the `1_NetColoc` notebooks
import, never the `2_CrossSpeciesBMI` ones. The only sensible read is that this was the single
conda environment (`cross_species_bmi`) the analyst used for both stages in sequence — NetColoc
first, then CrossSpeciesBMI in the same env, named after the later deliverable.

**One exception**: `1_NetColoc`'s `PLOT_NETCOLOC_SCATTER` process (Fig. 1g, added 2026-09-21) is
pure R/ggplot2 — this container is Python-only, so that one process reuses the shared
`gwas_supplement_plots` container instead (see
[`env/shared/README.md`](../shared/README.md)), not this one.

## Environment Files

### [`1_NetColoc/netcoloc.yml`](./1_NetColoc/netcoloc.yml)

Container for the `1_NetColoc`/`2_CrossSpeciesBMI` compute chain — network propagation z-scores
(`netcoloc`), community-detection API access (`cdapsutil`), NDEx interactome access (`ndex2`), and
supporting analysis packages.

| Tool | Version | Purpose |
|------|---------|---------|
| python | 3.9.13 | Runtime (`CrossSpeciesBMI`'s own pin) |
| netcoloc | 0.1.6.post1 | Network propagation z-scores (`netprop_zscore`, `netprop`, `network_colocalization`) |
| cdapsutil | 0.2.0a1 | Community-detection-as-a-service client (HiDeF, via the live CDAPS REST service) |
| ndex2 | 3.5.0 | Read/write access to NDEx-hosted networks (PCNet2.0 fetch) |
| mygene | 3.2.2 | Gene ID mapping, used in the `CrossSpeciesBMI` control analyses |
| networkx | 3.0 | Graph library — see version-conflict note below |
| numpy / pandas / scipy / statsmodels | 1.24.2 / 1.5.3 / 1.10.1 / 0.13.5 | Standard scientific stack, all `CrossSpeciesBMI`'s own exact pins |
| python-igraph | 0.10.4 | Graph library used by `ddot`/`cdapsutil` internals |
| tulip-python | 5.7.0 | Graph visualization library, a real (not incidental) `ddot` dependency |
| ddot (git) | tag `v1.0` | Data-Driven Ontology Toolkit — `michaelkyu/ddot`, installed via `--conda-run-command`, not pip |
| ndex-dev (pip, `--no-deps`) | latest | Older NDEx client `ddot` itself needs (`ndex.networkn.NdexGraph`) — see caveats below |

Validated with `pixi lock --dry-run` against `conda-forge` for `linux-64` — resolves cleanly.
Confirmed in the built container: every package above imports successfully, including the full
`ddot` → `ndex-dev` → `ndex.networkn` chain.

**Two build issues in this environment, not visible from the package names alone:**

1. **`ddot==1.0` isn't PyPI's own "ddot" package** (an unrelated, different tool whose own
   releases only go up to `0.4.0`) — it's `michaelkyu/ddot` at git tag `v1.0`, whose `setup.py`
   reports `version='1.0'`, matching what `CrossSpeciesBMI`'s `pip freeze` captured from a
   git-installed copy. Installed via `git clone` + `pip install --no-deps`, same pattern as
   PolyFun in `env/04_gwas/3_finemap_susie/`.

2. **Both `ddot` and `ndex-dev`'s own `setup.py` pin the ancient `networkx==1.11`** — installing
   either without `--no-deps` silently downgrades the environment's `networkx` and breaks it
   outright (`ImportError: cannot import name 'gcd' from 'fractions'` — `fractions.gcd` was
   removed in Python 3.9, `networkx==1.11` predates the `math.gcd` migration). `CrossSpeciesBMI`'s
   real historical environment explicitly overrode this to `networkx==3.0`; this env file does the
   same, pinning `networkx==3.0` directly and installing both `ddot` and `ndex-dev` with
   `--no-deps` so neither can reintroduce the stale pin.

**`ndex-dev`'s only published wheel is tagged `py2-none-any`** (Python 2 only) — its `setup.py`
classifiers list only `Python :: 2.7`, and its whole release history (`3.0.11.2` through
`3.0.11.41`) never had a proper Python 3 wheel cut. The module `ddot` actually imports
(`ndex/networkn.py`) and its full import chain (`ndex/__init__.py`, `ndex/create_aspect.py`,
`ndex/client.py`) parse cleanly under Python 3, and `networkn.py` has an explicit
`try: basestring; except: basestring = str` compatibility shim plus `from six import
string_types` — deliberate Python 2/3 compatibility, not accidental. The `py2` wheel tag is stale
packaging metadata (built once under Python 2, never re-cut), not a reflection of the source.

**`matplotlib-venn` deliberately dropped from this pin set** (present in `CrossSpeciesBMI`'s real
environment at `0.11.9`) — it has no wheel on PyPI, sdist-only, which breaks `pixi`'s `--dry-run`
precheck (needs an installed environment to build from source, not just index metadata). Not
needed by anything in `1_NetColoc`/`2_CrossSpeciesBMI`'s compute chain — only by the deferred
plotting notebooks (`venn_network_gene_overlap.ipynb`). Add back when that piece is implemented.

### [`5_GO_semantic_clustering/go_semantic_plots.yml`](./5_GO_semantic_clustering/go_semantic_plots.yml)

One container for all 8 processes in `5_GO_semantic_clustering`
(`SEMANTIC_CLUSTERING`, `PLOT_GO_CLUSTERS`, `OCD_NA_DETAIL`, `PLOT_SEMANTIC_SIMILARITY`,
`PLOT_FOLD_ENRICHMENT`, `PLOT_GENE_ENRICHMENT_GRIDS`, `PLOT_GENE_ENRICHMENT_GRIDS_V2`,
`PLOT_GENES_FOR_MULTIPLE_TERMS`), also reused by `4_Plotting` (`PLOT_VENN_HIERARCHY`, `PLOT_FIG2A`)
and `07_additional_plots_and_analyses/4_power_analysis` — they all share the same package list:
`r-base`, `r-tidyverse`, `r-glue`, `r-httr2`, `r-jsonlite`, `r-googlesheets4`, `r-cowplot`,
`r-tidytext`, `r-scales`, `r-ape`, `r-cluster`, `bioconductor-gosemsim`, `bioconductor-org.hs.eg.db`,
`bioconductor-go.db`, `bioconductor-annotationdbi` — all conda-forge/bioconda, no source builds
needed.

**`r-cluster` added and rebuilt 2026-10-03** for `01_semantic_clustering.R`'s updated
cutoff-sensitivity scan (`cluster::silhouette()`) — not in the originally frozen image. Rebuilt
with the same command below; the new tag is updated everywhere it's referenced (this file and
every process in `5_GO_semantic_clustering`, `4_Plotting`, and `4_power_analysis` that reuses this
container — all three were on the same tag before and after, so no process lost the packages it
already had).

Precheck (`pixi lock --dry-run` against the spec, no install) confirmed the environment resolves
before building. Built with the Wave CLI (v1.8.2):

```bash
wave --conda-file env/05_network/5_GO_semantic_clustering/go_semantic_plots.yml --platform linux/amd64 --freeze --await -o json
```

Originally frozen to `community.wave.seqera.io/library/go_semantic_plots:2df904e3debeab8d`. Not run
against real data end-to-end from inside the built container in that original pass — see
`archive/conversion_notes/05_network.md` for why (no R/tidyverse available locally to drive that
validation; only the DSL2 wiring was stub-tested). This stage was originally built (and its
container frozen) under `07_additional_plots_and_analyses`; moved here 2026-09-15.

**Rebuilt 2026-10-03** (same command, updated spec with `r-cluster` added) to
`community.wave.seqera.io/library/go_semantic_plots:69629d5c597a400b` — the current tag everywhere
this container is referenced. Not independently re-validated end-to-end against real data (the
underlying script's own logic was verified separately, byte-for-byte, against the collaborator's
real output — see `code/05_network/5_GO_semantic_clustering/README.md`'s Notes); this build only
adds one package on top of an already-working image.

### [`4_Plotting/network_overlap.yml`](./4_Plotting/network_overlap.yml)

Container originally built for `PLOT_NETWORK_OVERLAP` (Fig. 2a-b, added 2026-09-22) and
`2_CrossSpeciesBMI`'s `PLOT_CONSERVED_NETWORK_VENN` (Fig. 1f, added 2026-09-22 — both were simple
`eulerr` Venn/plotting jobs, no new packages needed): `r-base`, `r-ggplot2`, `r-dplyr`, `r-tidyr`,
`r-patchwork`, `r-eulerr`. `eulerr` wasn't pinned anywhere else in this repo before this, so this
got a dedicated build rather than reusing an unrelated container.

**2026-10-01: `4_Plotting`'s Panel B no longer uses this container** — its Euler diagram was
replaced by area-proportional Venn diagrams with no `eulerr` dependency (`PLOT_VENN_HIERARCHY`,
reuses `go_semantic_plots` instead; see `code/05_network/4_Plotting/README.md`'s Notes). This
container is still live for `2_CrossSpeciesBMI`'s `PLOT_CONSERVED_NETWORK_VENN` (Fig. 1f), which is
unaffected.

Precheck (`pixi init --import` + `pixi lock --dry-run`, no install) confirmed the spec resolves
cleanly before building. Built with the Wave CLI (v1.8.2):

```bash
wave --conda-file env/05_network/4_Plotting/network_overlap.yml --platform linux/amd64 --freeze --await -o json
```

Frozen to `community.wave.seqera.io/library/network_overlap:282efeb679f77cf4`. Validated end-to-end
(when still used by `4_Plotting`): pulled the built image and ran the archived
`plot_network_overlap.R` (now `archive/fig2b_euler_diagram_superseded/`) inside it against the real
deposited gene lists — produced the same output as the local run, matching the then-published
Fig. 2a-b numbers exactly.

## Notes

- All containers are built for `linux-64` platform.
- Channels: `conda-forge`.
