# Conversion notes — 05_network

Decisions and discrepancies noted while converting this stage's original scripts/notebooks into
the Nextflow pipelines under `code/05_network/`. Moved out of the `.nf`/`process_definitions.nf`
header comments and `bin/*` script docstrings so those files document only the science/pipeline
itself.

## code/05_network/1_NetColoc/process_definitions.nf

```
The notebooks recompute the SAME 4 gene sets' z-scores redundantly across separate,
hand-rerun kernel sessions (dogCD's z-scores alone get computed 3x across the 3 cross-species
pairings in 1.1, plus once more in 1.2) — an artifact of the notebook's own interactive,
cell-by-cell structure, not a real requirement. This pipeline computes each of the 4 distinct
gene sets (OCD, dogCD/CCD, DEP, SCH) exactly once; downstream cross-species combination steps
(2_CrossSpeciesBMI) reuse these same z-score outputs by reference.

Two of these 4 gene-set computations are FULLY reproducible and validated by real execution
against NetColoc's own source (not assumed):
  - w_prime/individual_heats_matrix (get_normalized_adjacency_matrix/get_individual_heats_matrix):
    pure deterministic linear algebra, no randomness anywhere. Recomputed live here rather than
    requiring the deposited 2.8GB precomputed matrices the original README offers as a shortcut
    — genuinely cheap (a dense O(N^3) inversion, a few minutes on one node with decent BLAS).
  - calculate_heat_zscores: has a `random_seed=1` DEFAULT parameter, called via
    `np.random.seed(random_seed)` as the function's very first line. Neither notebook ever
    overrides it, so every historical z-score run used this same fixed seed — reproduced
    explicitly here (not left to depend on netcoloc's own future default).

New container required: both NetColoc and CrossSpeciesBMI (github.com/ucsd-ccbb/NetColoc,
github.com/sarah-n-wright/CrossSpeciesBMI) are external Python tools. Built from
CrossSpeciesBMI's own real, fully-pinned environment.yml — see env/05_network/README.md for the
full build story, including two real build bugs found and fixed (ddot's true git source; a
networkx version conflict between ddot/ndex-dev's stale setup.py pins and the real environment).

[...]

One task per human trait paired against dogCD/CCD (OCD, DEP, SCH) — see bin/combine_zscores.py's
header for the source cell and the "r"/"h"/"hr" column-naming story. Validated by exact
recomputation of all 3 deposited *_zcomb_z12* files in data/05_network/ (numerically identical
to floating-point text-formatting precision — R's write.table truncates to ~15 significant
digits where pandas' default doesn't, so a byte-for-byte diff isn't expected, but np.allclose
confirms max abs diff ~1e-14 across all 3 pairings).
```

## code/05_network/2_CrossSpeciesBMI/2_CrossSpeciesBMI.nf

```
dogCD/CCD-OCD conserved network, colocalization size vs. a permuted null (Figure 1g), the
control-comparisons bar plot, the dogCD/CCD-OCD systems map hierarchy, and per-pairing MGD
significant-community derivation. See process_definitions.nf's header for the important caveat
about this notebook's missing local helper modules (now resolved — the real files are committed
alongside the notebooks) and the `final_annotations` gap.

[...]

All 7 systems-map pairings, each with a real, correctly-formatted deposited MGD enrichment
result — see process_definitions.nf's header for why this is generalized across all of them
rather than just the OCD-CCD pairing 2.2's own cells literally show, matching this project's
established "genericize the worked example" practice. OCD_ONLY's own deposit originally had
unrelated LDSC annotation content, not MGD enrichment results — the collaborator re-deposited it
(main branch commit 34ccf40, as fixedOCD_ONLY_hierarchy_full_MGD_enrichment_results.csv, in a
different export format than its siblings); reformat_ocd_only_mgd_results.py normalizes it to
match, run once to produce the OCD_ONLY_hierarchy_full_MGD_enrichment_results.tsv referenced here.
```

## code/05_network/2_CrossSpeciesBMI/process_definitions.nf

```
This notebook (and 2.2/2.3) import 3 local helper modules — analysis_functions.py,
plotting_functions.py, updated_netcoloc_functions.py — via `sys.path.append(cwd); from X import
*`, now committed alongside the notebooks in this directory (the collaborator's own local
copies). Every function used here (calculate_network_overlap, calculate_expected_overlap,
plot_permutation_histogram, get_p_from_permutation_results) is transcribed directly from these
real files — confirmed byte-for-byte against them after they landed (this pipeline was
originally built against github.com/sarah-n-wright/CrossSpeciesBMI's public main-branch versions
as a stand-in; plotting_functions.py/updated_netcoloc_functions.py turned out identical, and
analysis_functions.py differed by exactly the one already-known, already-matched line — load_pcnet()'s
interactome UUID, overridden locally to PCNet2.0).

[...]

Observed network overlap size is fully deterministic; the permuted null distribution is NOT
(unseeded random.shuffle in the source) — a fresh permutation is generated each run, matching
this project's established treatment of every other non-reproducible permutation step.

[...]

controlanalyses_260521.csv is an already-deposited, frozen historical result (the 12-
combination control loop that produces it uses an unseeded permutation and is not
reproducible; the notebook's own cell 67 discards its freshly computed result and reloads
this same file before plotting, so this process does the same rather than re-deriving it).

[...]

Reloads the already-deposited, frozen MGD enrichment result rather than recomputing it — see
this file's header (2.2's own computation is commented out in favor of this same reload; 2.3's
isn't, but its own output filename is already deposited too, so the same "compute once, reuse"
pattern is applied consistently rather than re-running community_term_enrichment, which would
additionally need the full MGI/MPO/ontology machinery this pipeline doesn't otherwise build).
2 of the 6 deposited files (SCH_CCD, DEP_CCD) use a comma decimal separator instead of a
period — coerced explicitly in bin/mgd_significant_communities.py, see its own header.
```

## code/05_network/4_Plotting/process_definitions.nf

```
Output filenames are passed in explicitly (not derived from trait_id) because the original
notebook's own naming is inconsistent across traits (e.g. "GO_enrichment_CCDOCD_OCD_..." vs
"GO_enrichment_CCDDEP-DEP_..." — underscore vs hyphen) — reproduced literally rather than
normalized. Validated by exact byte-for-byte reproduction of all 3 deposited
*_GO_difference_volcano_all_terms.txt and *_run1_only_terms.txt files in
data/05_network/plotting/.
```

## code/05_network/1_NetColoc/bin/fetch_pcnet2.py

```
Source: 04_gwas.../1_NetColoc/1.1_OCD_dogCD_NetColoc_analysis_260521.ipynb, cell 10.

Literal reproduction of the notebook's own interactome-fetch cell, not NetColoc's own
netprop_zscore() wrapper function — the notebook fetches the interactome and calls
netprop_zscore.calculate_heat_zscores() directly, never the wrapper, so the wrapper's own
internal cleanup step ("remove 'None' node") is NOT what actually ran historically. The
notebook's own literal cleanup is "remove self-loop edges" (`G_int.remove_edges_from(
nx.selfloop_edges(G_int))`) — reproduced here instead.

interactome_uuid defaults to PCNet2.0 (d73d6357-e87b-11ee-9621-005056ae23aa, 19,267 nodes /
3,852,119 edges per the notebook's own comment) — the notebook also has PCNet (18,820/2,693,109)
and PCNet2.2 (18,558/3,323,928) commented out as alternatives it tried; PCNet2.0 is what was
actually left active and used.

The raw CX export for this network is >1GB (confirmed via a direct curl check, not assumed) —
much larger than the summarized "19,267 nodes" figure suggests, since CX format repeats node/edge
attributes verbosely. Real testing surfaced intermittent `IncompleteRead`/`ChunkedEncodingError`
failures partway through the streamed download (a genuine, reproducible network-reliability issue
for a transfer this large, not a bug in this script) — retried here with backoff rather than
relying on Nextflow's own task-level retry, which would re-run the entire task just to retry one
HTTP call.
```

## code/05_network/1_NetColoc/bin/build_heats_matrix.py

```
Source: 1.1_OCD_dogCD_NetColoc_analysis_260521.ipynb, cell 13 (the commented-out recompute code —
the notebook itself loaded CrossSpeciesBMI's already-deposited matrices instead: "We loaded the
ones created in CrossSpecies instead"). No literal historical run of this recompute path exists,
but the code is the original author's own, just never executed in this exact notebook run — it was
executed once, earlier, to produce the deposited .npy files this notebook then reused.

Deliberately recomputed live here rather than requiring the deposited 2.8GB matrices: both
get_normalized_adjacency_matrix/get_individual_heats_matrix are pure deterministic linear algebra
(no randomness anywhere) — checked directly against NetColoc's own source. The only cost is a
dense O(N^3) matrix inversion (~7x10^12 FLOPs for N=19,267), a few minutes on one HPC node with a
decent BLAS backend, ~10-15GB peak RAM — matching the notebook's own comment ("this step takes a
few minutes, more for denser interactomes"). Recomputing avoids depending on an external 2.8GB
deposit at all.
```

## code/05_network/1_NetColoc/bin/combine_zscores.py

```
Source: 1.1_OCD_dogCD_NetColoc_analysis_260521.ipynb, cell 24 ("in R:" block) — full_join(ccd, ocd,
by="gene") then zcomb = D1_z * D2_z, columns renamed to NPS_r/NPS_h/NPS_hr (per 2_CrossSpeciesBMI's
own cell-10 comment: "h = human, r = dog, hr = dog-human"; "r" is a historical holdover from an
earlier version of this pipeline that compared against rat, kept for column-naming consistency
across all 3 deposited *_zcomb_z12* files even though every non-human trait here is the dog CCD set).

No equivalent code cell exists in the repo for the DEP-CCD/SCH-CCD pairings (only OCD-CCD's is
shown, cell 24), but DEP_CCD_Oct25_zcomb_z12.txt/SCH_CCD_zcomb_z12.txt are already deposited in
data/05_network/ following the exact same 3-column mechanic (confirmed: their NPS_r column is
byte-identical to ocd_ccd_zcomb_z12_251022.txt's, since CCD is always the "r" side) — generalized
here by argument rather than hand-copied per pair, same practice as compute_zscores.py.

R's full_join keeps every gene present in either input (outer join), with NA where one side is
missing; pandas' outer join + multiplication does the same (NaN propagates through '*' exactly like
R's NA). Row order follows the "r" (dog/CCD) input first, matching full_join(ccd, ocd, ...)'s own
argument order in the source cell.
```

## code/05_network/1_NetColoc/bin/compute_zscores.py

```
Source: 1.1_OCD_dogCD_NetColoc_analysis_260521.ipynb, cells 6/7 (seed gene loading), 15 (interactome
intersection), 17/20 (the D1/D2 z-score calls — identical shape, genericized here across all 7
traits by argument rather than by hand-editing per pair, same "genericize the worked example"
practice as every other stage in this project).

Seed genes are loaded via pandas (a single-column CSV with a header named "gene"), NOT the plain
whitespace-split text file netcoloc's own top-level netprop_zscore() wrapper expects — the notebook
never calls that wrapper; it calls calculate_heat_zscores() directly and does its own seed-gene
loading and interactome intersection first (cell 15's intersection is redundant with
calculate_heat_zscores()'s own internal np.intersect1d call, but reproduced anyway for fidelity).

alpha=0.5 and random_seed=1 are never written explicitly in the notebook's own call (cells 17/20
only override num_reps=1000 and minimum_bin_size=100) — both are the function's own defaults.
Passed explicitly here as script arguments (still defaulting to the same values) so the exact
values used don't silently depend on whatever netcoloc's own defaults happen to be in a future
release.

Output matches the deposited historical z-score CSVs exactly: a single 'z' column (from the
notebook's own `pd.DataFrame({'z': z_D1})`), written with the default pandas index (gene names).
```

## code/05_network/2_CrossSpeciesBMI/bin/conserved_network.py

```
Source: 2.1_OCD_dogCD_Network_Colocalization_260521.ipynb, cell 26. Reproduces the notebook's
literal manual threshold (z_human > 1.5 & z_rat > 1.5 & z_human*z_rat > 3) directly, rather than
calling CrossSpeciesBMI's own `calculate_network_overlap` helper — mathematically equivalent here
since both individual thresholds already force z_human/z_rat positive, making that helper's extra
`*(z1>0)*(z2>0)` masking a no-op for this specific threshold combination (confirmed by reading its
source: sarah-n-wright/CrossSpeciesBMI updated_netcoloc_functions.py).

The notebook's own custom helper modules (analysis_functions.py, plotting_functions.py,
updated_netcoloc_functions.py) are committed alongside the notebooks in this directory (the
collaborator's real local copies). This script's logic was confirmed byte-for-byte against
updated_netcoloc_functions.py's real calculate_network_overlap/calculate_expected_overlap after
they landed.
```

## code/05_network/2_CrossSpeciesBMI/bin/control_analysis_plot.py

```
Source: 2.1_OCD_dogCD_Network_Colocalization_260521.ipynb, cell 70. Consumes the already-deposited
`controlanalyses_260521.csv` (data/05_network/) — the 12-combination control loop that produces it
(cell 65) calls `calculate_expected_overlap` with an unseeded permutation (num_reps=1000, no
overlap_control), so it is not reproducible; the notebook's own cell 67 discards its freshly
computed `control_results` and reloads this same frozen CSV before plotting, so this script does
the same rather than re-deriving it.

DISCREPANCY: cell 67 reads this file with `sep=","`, but the actual deposited
controlanalyses_260521.csv is semicolon-delimited (confirmed by inspection) — read with `sep=";"`
here to make it parse at all; flagged, not silently worked around.
```

## code/05_network/2_CrossSpeciesBMI/bin/merge_hierarchy.py

```
Source: 2.2_OCD_dogCD_Systems_Map_260521.ipynb, cells 12-17. Literal reproduction: `hier_df_genes`
(CD_MemberList, represents, indexed by community name) from the frozen hierarchy data, joined
(inner) against `frac_d1_seeds`/`frac_d2_seeds` from the systems map, matching cell 16's
`pd.merge(..., left_index=True, right_index=True, how='inner')` and cell 17's dtype coercions.
```

## code/05_network/2_CrossSpeciesBMI/bin/mgd_significant_communities.py

```
Source: 2.2_OCD_dogCD_Systems_Map_260521.ipynb / 2.3_OCD_ONLY_Systems_Map_260525.ipynb, cells 36-42
(2.2's own MGD enrichment computation is commented out in favor of reloading the deposited result;
2.3's isn't commented out, but its own output filename IS already deposited — same "compute once,
reuse" pattern applied consistently here, so this script reloads rather than re-running
`community_term_enrichment`, which additionally needs the full MGI/MPO/ontology machinery this
pipeline doesn't otherwise build — see this directory's README for why that's out of scope).

DELIBERATELY SKIPS the `.merge(final_annotations, ...)` step in both notebooks' cell 37:
`final_annotations` is never defined anywhere in either notebook (checked every cell in both) — it
looks like a lost cell from whatever produced the deposited notebooks. Investigated whether it's
actually needed: the "description" column both notebooks report via `value_counts("description")`
is already present in the raw MGD enrichment file itself (one row per (community, MP-term) pair,
each carrying that MP term's own description, e.g. "abnormal behavior") — confirmed by inspecting
the file directly, not assumed. Nothing downstream of the merge (cells 37, 39-42) touches any column
that could only come from `final_annotations`, so this script uses the base MGD file's own
"description" column throughout and omits the merge rather than guessing at its content.

DISCREPANCY: 2 of the 6 deposited MGD files (SCH_CCD, DEP_CCD) use a comma as the decimal separator
in their numeric columns (e.g. "2,41E-116"), unlike the other 4 (period-decimal) — confirmed by
inspection, not assumed. Read naively, this silently loads OR_p/OR_CI_lower etc. as strings, which
both crashes fdrcorrection() AND makes sort_values(by="OR_CI_lower") do a lexicographic string sort
instead of a numeric one (a silent correctness bug, not just a crash). Coerced explicitly below,
flagged rather than silently normalized.
```

## code/05_network/2_CrossSpeciesBMI/bin/netcoloc_size_permutation.py

```
Source: 2.1_OCD_dogCD_Network_Colocalization_260521.ipynb, cells 28/30 — calls
`calculate_expected_overlap(z_rat, z_human, z_score_threshold=3, z1_threshold=1.5, z2_threshold=1.5,
num_reps=10000, overlap_control="bin", seed1=dog_seeds, seed2=human_seeds)` then
`plot_permutation_histogram(permuted, observed, ...)`.

`calculate_expected_overlap` is a CrossSpeciesBMI-specific reimplementation (NOT netcoloc's own
`network_colocalization` module, which this notebook doesn't even import) — its logic here was
confirmed byte-for-byte against the collaborator's real updated_netcoloc_functions.py/
plotting_functions.py, now committed alongside the notebooks in this directory.

REPRODUCIBILITY: the "observed" network overlap size is fully deterministic (computed once from the
real z-scores, no randomness) and IS reproduced exactly. The permuted null distribution is NOT
reproducible — the source shuffles with `random.shuffle()` and never seeds it — so this script
generates a genuine, fresh permutation each run (matching this project's "frozen data architecture"
treatment of every other non-reproducible permutation/HiDeF step: executable and real today, just
not byte-identical to the historical run). The resulting p-value will therefore vary slightly
run-to-run around the same true value, not because of a bug.
```

## code/05_network/4_Plotting/bin/go_term_comparison.py

```
Source: 05_network/4_Plotting/GO-term_comparison.ipynb, cells 1/2/3 — identical logic across the 3
gene-set pairs (OCD, DEP, SCH), each comparing a single-species-only hierarchy GO enrichment
("run1") against the corresponding dogCD/CCD cross-species one ("run2").

xlim/ylim: cell 1 (OCD) computes the plot limits dynamically from the data; cells 2/3 (DEP, SCH)
hardcode them instead (`plt.xlim(-20,70)` / `plt.ylim(0,140)`, with the author's own comment "change
to match schizophrenia" — a deliberate visual-comparison choice, not an oversight) — reproduced
literally via optional --xlim/--ylim rather than "fixed" to be consistent across all 3 traits.
```

(Note: this entry duplicates the `code/05_network/4_Plotting/process_definitions.nf` heading above
for the GO_TERM_COMPARISON process comment — both the process comment and this script's own
docstring carried overlapping narration and are archived from their respective files.)

## env/05_network/README.md

```
**Two real build issues found and fixed, not just assumed away:**

1. **`ddot==1.0` isn't PyPI's own "ddot" package** (an unrelated, different tool whose own
   releases only go up to `0.4.0`) — it's `michaelkyu/ddot` at git tag `v1.0`, whose `setup.py`
   reports `version='1.0'`, matching exactly what `CrossSpeciesBMI`'s `pip freeze` captured from a
   git-installed copy. Installed via `git clone` + `pip install --no-deps`, same pattern as
   PolyFun in `env/04_gwas/3_finemap_susie/`.

2. **Both `ddot` and `ndex-dev`'s own `setup.py` pin the ancient `networkx==1.11`** — installing
   either without `--no-deps` silently downgrades the environment's `networkx` and breaks it
   outright (`ImportError: cannot import name 'gcd' from 'fractions'` — `fractions.gcd` was
   removed in Python 3.9, `networkx==1.11` predates the `math.gcd` migration). `CrossSpeciesBMI`'s
   real historical environment explicitly overrode this to `networkx==3.0`; this env file does the
   same, pinning `networkx==3.0` directly and installing both `ddot` and `ndex-dev` with
   `--no-deps` so neither can reintroduce the stale pin. Confirmed via a real, non-stub container
   run — the first build attempt (without this fix) failed with exactly that `ImportError`.

**`ndex-dev`'s only published wheel is tagged `py2-none-any`** (Python 2 only) — its `setup.py`
classifiers list only `Python :: 2.7`, and its whole release history (`3.0.11.2` through
`3.0.11.41`) appears to have never had a proper Python 3 wheel cut. Checked before trusting a
source build under Python 3.9 would even work: the actual module `ddot` imports
(`ndex/networkn.py`) and its full import chain (`ndex/__init__.py`, `ndex/create_aspect.py`,
`ndex/client.py`) all parse cleanly under Python 3, and `networkn.py` has an explicit
`try: basestring; except: basestring = str` compatibility shim plus `from six import
string_types` — deliberate Python 2/3 compatibility, not accidental. The `py2` wheel tag looks like
stale packaging metadata (built once under Python 2, never re-cut), not a reflection of the
source. Confirmed empirically, not just by static reading: the container builds and imports the
full chain successfully.
```

## 5_GO_semantic_clustering

Decisions and discrepancies noted while converting a collaborator's GO semantic-clustering /
fold-enrichment / community-tileplot chain into the Nextflow pipeline under
`code/05_network/5_GO_semantic_clustering/`. Every fix described here is a Nextflow-level
input-staging decision, not a script-logic change — none of the 7 scripts' actual analysis logic
has ever been edited (see the 2026-09-15 reorg note directly below for the 2 narrow, direct-
instruction exceptions: a staging-directory rename and a shebang/execute-bit addition, neither of
which touches what the scripts compute). (Built under `07_additional_plots_and_analyses` on
2026-09-11; moved to `05_network` on 2026-09-15 per KatarinaTe's request to keep it alongside the
rest of the network analyses rather than split off on its own.)

**bin/ reorganisation and further corrections (2026-09-15), per direct instruction.** The general
"at the end of this refactor, only Nextflow-related files needed to produce the paper's results
stay outside archive/" policy applied here first: the 7 scripts moved from the collaborator's
original working directory into this stage's own `bin/` (matching every other stage's convention —
no longer passed in as `path(script)` process inputs/pipeline params), the 2 frozen inputs
(`cache_gpt_labels.top_community.cut60.txt`, `goa_human.gaf.gz`) moved up to sit directly under
`data/05_network/5_GO_semantic_clustering/`, and everything else that had been in that original
working directory — split by which of two unrelated analyses each file actually belonged to — moved
to `archive/go_semantic_clustering_originals/` (pre-generated figure outputs, `.Rproj`/review files)
and `archive/gwas_catalog_enrichment_originals/` (the separate GWAS-catalog "12 of 27" claim's
scripts/cache). Per further direct instruction: the 7 scripts were each given an execute bit and (4
of them; 3 already had one) a prepended `#!/usr/bin/env Rscript` shebang, so they run bare via
Nextflow's automatic `bin/`-on-`PATH`, exactly like every other stage's `bin/` scripts, rather than
via an explicit `Rscript "<path>"` call; and the `"../Figure4/..."` staging directory two of the
scripts hardcode (`01_semantic_clustering.R`, `plot_fold_enrichment.R`) was renamed to
`"../go_enrichment_tsvs/..."` in both the scripts and `process_definitions.nf`'s matching staging —
"Figure4" tied the directory's name to a manuscript figure number that could change, when its actual
content (6 GO enrichment TSVs) doesn't depend on that numbering at all. These are the two narrow
exceptions to "content never touched" for this stage, both authorized directly rather than assumed.

That same pass also corrected a real, previously wrong attribution, found while renaming the
staging directory and re-checking which scripts actually use it: **`plot_fold_enrichment.R`, not
`plot_gene_enrichment_grids.R`, is the actual builder of `all_go_combined.rds`** — see the
`PLOT_GENE_ENRICHMENT_GRIDS`/`PLOT_FOLD_ENRICHMENT` entries below for the corrected reasoning; the
original investigation (and an intermediate, also-wrong attempt earlier the same day to model it as
a missing/frozen input) are left in place below rather than deleted, so the correction is visible in
context.

### General staging strategy (process_definitions.nf)

Every script here is a plain `Rscript` with no CLI arguments: settings and output filenames are
hardcoded constants, and file reads use paths relative to whatever directory the script happens to
be run from. Rather than rewrite any of this logic into a parameterized bin/ script, each Nextflow
process stages the exact inputs each script expects at the exact relative path it expects them, then
runs the unmodified script. Two scripts (`01_semantic_clustering.R`, `plot_fold_enrichment.R`) read
from `"../go_enrichment_tsvs/<file>.tsv"` (named `"../Figure4/..."` until the 2026-09-15 rename
above) relative to the script's own working directory — handled by running each of those 2 scripts
one directory below a `go_enrichment_tsvs/` directory populated with the 6 real GO enrichment TSVs
(see below), inside the task's own work directory.

### `../Figure4/...` broken path, later renamed `../go_enrichment_tsvs/...` (01_semantic_clustering.R,
### plot_fold_enrichment.R)

Both scripts read 6 GO-enrichment TSVs from a `../Figure4/` directory that does not exist anywhere
in this repo. Confirmed by `diff`: the real files are at `data/05_network/<same filename>.tsv` for
all 6 — `DEP_ONLY_hierarchy_full_GO_enrichment.tsv`, `DEP_CCD_hierachy_full_GO_enrichment.tsv`,
`OCD_ONLY_hierachy_full_GO_enrichment.tsv`, `Compulsive_hierachy_full_GO_enrichment_251023.tsv`,
`SCH_ONLY_hierachy_full_GO_enrichment.tsv`, `SCH_CCD_hierachy_full_GO_enrichment.tsv`. This
mismatch is already flagged in `REPRODUCIBILITY_AUDIT.md`'s 2026-09-11 update; the initial pass
fixed it as a Nextflow staging decision (a `Figure4/` directory built inside the task work dir), not
a script edit; the 2026-09-15 pass (see above) went further and renamed the directory itself, in
both `process_definitions.nf` and the 2 scripts, to `go_enrichment_tsvs/` — content-descriptive
rather than tied to a manuscript figure number.
(`plot_gene_enrichment_grids.R` was originally, incorrectly, included in this "3 scripts" — see its
own entry below for the correction: it never reads these TSVs at all.)

### `01_semantic_clustering.R`

**Unconditional `OPENAI_API_KEY` check.** The script's GPT-labeling section starts with
`if (Sys.getenv("OPENAI_API_KEY") == "") stop(...)`, evaluated *before* checking whether the
frozen `cache_gpt_labels.top_community.cut60.txt` already covers every cluster that needs a label.
Against a clean rerun of the frozen inputs, the anti-join against that cache is expected to leave 0
clusters needing a fresh GPT call — but the env-var check still fires regardless, so a real key is
never actually used in that case. `params.openai_api_key` defaults to a placeholder string
(`'unused-frozen-cache-only'`) passed as `OPENAI_API_KEY` purely to satisfy this check. If the
semantic-clustering step (recomputed fresh every run — its own `.rds` similarity-matrix cache is
NOT checked in) produces cluster IDs the frozen cache doesn't already have a label for, a real
OpenAI API key would be genuinely required for those specific new clusters; not fabricated here.

**Final `googlesheets4::sheet_write()` call.** The script's very last substantive action mirrors
its per-term output table to a Google Sheet tab, with no `gs4_deauth()` call anywhere in this
script (unlike the group-4 tileplot scripts, which do call it before their own, read-only sheet
access) — meaning this call needs a cached, interactive OAuth token that does not exist in a
non-interactive container. This call runs strictly after both of the script's real, needed TSV
outputs (`GO_clusters_labeled.cut60.txt`, `GO_with_without_dogCD.semantic_clustering.
top_community.cut60.txt`) are already written to disk. `SEMANTIC_CLUSTERING` therefore runs the
script with `set +e`, captures its exit code, and only fails the task if the two required TSVs are
missing — a non-zero exit is tolerated (and logged loudly to stderr) if they're both present. This
is the same category of live, non-reproducible Google-Sheet write already documented elsewhere in
this repo (Fig 1b's heritability sheet, the GWAS-catalog `gene_list` sheet) — flagged here as a
known, expected failure mode, not silently worked around by, e.g., deleting the call.

**A stray, unreproduced output file already sat in the collaborator's original working directory:**
`GO_with_without_dogCD.semantic_clustering.cut60.txt` (no `top_community` in the name). The
current script only ever writes `GO_with_without_dogCD.semantic_clustering.top_community.cut60.txt`
(the `"top_community"` substring is a literal, unconditional part of the `paste0(...)` that builds
`out_file_long` — not dependent on the `restrict_to_top` toggle). This means the non-`top_community`
file was produced by some other/earlier version of this script, not the one now checked in — flagged
at the time, not deleted (it and 2 similarly cut50/cut70-suffixed siblings were later moved to
`archive/elinors_analysis_superseded/` as draft/orphaned variants).

### `02_plot_GO_clusters.R`, `03_OCD_NA_detail.R`

Both read via `list.files(".", pattern = "^GO_with_without_dogCD\\.semantic_clustering\\.
top_community\\..*\\.txt$")` — a plain working-directory glob, no broken-path fix needed. Each
process stages `SEMANTIC_CLUSTERING`'s per-term output directly into its own working directory and
runs the script as-is.

### `plot_gene_enrichment_grids.R` vs. `.v2.R` vs. `plot_genes_for_multiple_terms.R`

All 3 share the same "gene x GO-term nesting tileplot" logic (the latter two even carry the
identical internal docstring header, `"go_nesting_tileplot.R"` — a stale filename comment, not
evidence they're duplicates of each other's current behaviour). They are NOT interchangeable:

- **CORRECTION (2026-09-15):** all three of `plot_gene_enrichment_grids.R`, `.v2.R`, and
  `plot_genes_for_multiple_terms.R` — including the non-`.v2` one, contrary to what this note
  originally said — only ever do `stopifnot(file.exists(cache_path)); readRDS(cache_path)`. None of
  them build `all_go_combined.rds` themselves, and none of them read the 6 raw `go_enrichment_tsvs/`
  (formerly `Figure4`) TSVs at all — confirmed by grepping all three files in full for every raw
  TSV filename and for `saveRDS`, with zero hits in any of them.
- **The actual builder is `plot_fold_enrichment.R`** — a 4th script, not one of these three. It
  computes `cache_is_stale` (true whenever `all_go_combined.rds` doesn't already exist, which it
  never does in a fresh Nextflow task directory) and, on that branch, reads the 6
  `go_enrichment_tsvs/` files and `saveRDS`s the result, as a side effect of building its own
  fold-enrichment barchart (Fig 3b). See its own process entry for the wiring.

This pipeline now wires that dependency correctly: `PLOT_FOLD_ENRICHMENT` emits `all_go_combined.rds`
as a process output, consumed by `PLOT_GENE_ENRICHMENT_GRIDS`, `PLOT_GENE_ENRICHMENT_GRIDS_V2`, and
`PLOT_GENES_FOR_MULTIPLE_TERMS`. None of those three need the `go_enrichment_tsvs/` staging at all —
only `PLOT_FOLD_ENRICHMENT` and `SEMANTIC_CLUSTERING` do.

### Which script maps to which checked-in deliverable

- `plot_gene_enrichment_grids.R` (v4/non-`.v2`): `communities <- c("C185","C187","C186","C197")`
  (4), output stem is **literally hardcoded** `paste0(plot_ds, "_all_communities_p4_maxts",
  max_termsize)` — i.e. the `"p4"` is a fixed string, NOT `-log10(pvalue_max)` (which is actually
  `1e-5` here, i.e. `5`, not `4`). This exactly matches the checked-in
  `OCD_CCD_all_communities_p4_maxts1000_nesting_tileplot.{pdf,_data.csv,_legend.info.txt}` and
  `OCD_CCD_community_gene_overlap.tsv` (4-community version) — reproduced as literally
  hardcoded, not "corrected" to `p5` to match the dynamic value. **Resolved with KatarinaTe
  (2026-09-11): confirmed a plain mislabeling**, not a meaningful version marker — no earlier
  revision of this script ever used `pvalue_max <- 1e-4` anywhere. Left as-is since it's
  original-author code.

  **No PNG is ever written by this script** (only `ggsave(out_pdf, ...)` — no matching `ggsave` for
  a `.png`), despite a `OCD_CCD_all_communities_p4_maxts1000_nesting_tileplot.png` already sitting
  in `archive/go_semantic_clustering_originals/` alongside the PDF. Flagged, not fabricated:
  `PLOT_GENE_ENRICHMENT_GRIDS` does not emit a PNG output.

- `plot_gene_enrichment_grids.v2.R`: `communities <- c("C184","C185","C187","C186","C189","C197",
  "C193")` (7) — matching the 7 communities actually present in `tileplots_pdfs_and_csvs/`
  (`C184,C185,C186,C187,C189,C193,C197`, confirmed by listing that directory). Its own output stem
  is `paste0(plot_ds, "_", cm, "_p", -log10(pvalue_max), "_maxts", max_termsize)` — computed
  dynamically, unlike v4's hardcoded `"p4"`. With this script's **current** `pvalue_max <- 1e-5`,
  `-log10(pvalue_max)` = `5`, so a fresh run produces filenames like
  `OCD_CCD_C185_p5_maxts1000_nesting_tileplot.pdf` — **not** `..._p4_...`, which is what's already
  checked into `tileplots_pdfs_and_csvs/`. **Resolved with KatarinaTe (2026-09-11): not a bug.**
  `p5` (the script's current, correct config) is the pipeline's real output; the checked-in `p4`
  files are an earlier snapshot that predates a manual Cytoscape-amendment step the `p5` output is
  meant to feed — the same category of terminal manual polish already documented elsewhere in this
  repo (Fig 2b/3c's Cytoscape+Illustrator steps). `PLOT_GENE_ENRICHMENT_GRIDS_V2`'s outputs are
  declared via a `OCD_CCD_*_nesting_tileplot.{pdf,_data.csv}` glob (not hardcoded to `p4` or `p5`),
  so this pipeline correctly reports whatever the script's current config computes (`p5`).

- `plot_genes_for_multiple_terms.R`: `communities <- c("C184","C185","C186","C193","C197")` (5,
  a strict subset of v2.R's 7), `p_cutoff <- 1e-4` (v2.R: `0.01`), `max_termsize <- 500` (v2.R:
  `1000`), and its own output stem, `paste0(plot_ds, "_", cm, "_p", -log10(pvalue_max))`, has
  **no** `"_maxts"` component at all — a structurally different filename shape from both other
  scripts. This script is NOT referenced anywhere in `data/Manuscript_draft/
  REVIEW_AND_SUGGESTIONS.md`'s figure-by-figure source table (which lists only
  `plot_gene_enrichment_grids.R`/`.v2.R` for Fig 3d/e) or in `REPRODUCIBILITY_AUDIT.md`. Its
  relationship to the other two — an earlier draft, an alternate analysis, or something else —
  could not be determined from the repo alone. Wired into this pipeline anyway (named explicitly
  in scope for this conversion), publishing to its own `genes_for_multiple_terms/` output
  directory, but flagged as an artifact whose provenance/purpose relative to the manuscript's
  actual figures is unresolved.

### `OCD_CCD_community_gene_overlap.tsv` name collision

All 3 scripts independently write `paste0(plot_ds, "_community_gene_overlap.tsv")` —
`"OCD_CCD_community_gene_overlap.tsv"` — but each over a *different* community subset (4, 7, and 5
communities respectively). Same filename, 3 different contents. Not fixable without editing the
scripts (which is off-limits), so each process's overlap TSV is published to its own output
subdirectory (`gene_enrichment_grids/`, `tileplots_pdfs_and_csvs/`, `genes_for_multiple_terms/`)
rather than a shared one, to avoid one silently overwriting another.

### Seed-gene Google Sheet (all 3 group-4 scripts)

`seed_sheet <- "<Google Sheet link removed>
edit?usp=sharing"`, tab `"gene list w gwas catalog overlap"`, read via `googlesheets4::read_sheet()`
after `googlesheets4::gs4_deauth()` — i.e. an anonymous, public read (no login), IF the sheet is
actually shared "anyone with the link can view". Checked via a live fetch of the sheet URL: the
result was ambiguous — a sign-in interstitial was reported, but the response also surfaced real
tab names (including `"gene list w gwas catalog overlap"`, the exact tab these scripts read) and
partial cell content. This could not be resolved to a definite yes/no without an actual R session
driving `googlesheets4::gs4_deauth() + read_sheet()` against it, which this environment cannot run.
Per this project's established practice for exactly this situation (Fig 1b's heritability sheet,
the GWAS-catalog `gene_list` sheet): **not fabricated** — the 3 processes that need this sheet
(`PLOT_GENE_ENRICHMENT_GRIDS`, `PLOT_GENE_ENRICHMENT_GRIDS_V2`, `PLOT_GENES_FOR_MULTIPLE_TERMS`)
run the unmodified `gs4_deauth()`/`read_sheet()` call as written; if the sheet turns out to require
authentication, this becomes a genuine, undocumented external blocker for all 3, same category as
this repo's other frozen/manual Google-Sheet inputs, requiring either sheet-sharing settings to be
fixed by its owner or the sheet's content to be frozen into a repo-tracked file the same way
`cache_gpt_labels.top_community.cut60.txt` already is for the GPT labels.

### `goa_human.gaf.gz`

All 3 group-4 scripts fall back to downloading this file
(`http://current.geneontology.org/annotations/goa_human.gaf.gz`) if it isn't already present in
the working directory. Since a copy is already checked into this repo (frozen, now at
`data/05_network/5_GO_semantic_clustering/goa_human.gaf.gz`), each process stages that frozen copy
directly instead, so no live download is attempted or needed.

### Redundant `all_go_combined.rds`-equivalent work not otherwise de-duplicated

`PLOT_GENE_ENRICHMENT_GRIDS` builds this cache once; `PLOT_GENE_ENRICHMENT_GRIDS_V2` and
`PLOT_GENES_FOR_MULTIPLE_TERMS` each still redo their own downstream filter/derive/community work
independently from it (each script's own logic, unmodified) — no attempt was made to further
de-duplicate or refactor across the 3 scripts, matching this project's stated preference not to
restructure for a marginal gain when it would mean touching read-only original-author code.

### Container (`env/05_network/5_GO_semantic_clustering/go_semantic_plots.yml`)

One container built for all 7 processes: `r-base`, `r-tidyverse`, `r-glue`, `r-httr2`,
`r-jsonlite`, `r-googlesheets4`, `r-cowplot`, `r-tidytext`, `r-scales`, `r-ape`,
`bioconductor-gosemsim`, `bioconductor-org.hs.eg.db`, `bioconductor-go.db`,
`bioconductor-annotationdbi` — all conda-forge/bioconda, no source builds. Precheck (`pixi lock
--dry-run`) confirmed the spec resolves; built for real with the Wave CLI (v1.8.2, checksum
verified) and frozen to
`community.wave.seqera.io/library/go_semantic_plots:2df904e3debeab8d`. Not run against real data
end-to-end from inside the built container in this pass (no R/tidyverse available locally to drive
that validation) — only the DSL2 wiring was stub-tested (`nextflow run ... -stub`).
