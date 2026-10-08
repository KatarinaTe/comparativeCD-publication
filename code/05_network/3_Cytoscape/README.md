# 05. Network — 3. Cytoscape

**Manual, terminal step — no Nextflow pipeline.** Network/systems-map visualization
(Cytoscape v3.10.2), producing the hand-laid-out figures (Fig. 2b, 3a, etc.) that Illustrator
then polishes. Not re-runnable by design — see `REPRODUCIBILITY_AUDIT.md`'s "Human-in-the-loop
checkpoints" for why this is an accepted, documented limitation rather than a gap.

## Files

| File | Description |
|---|---|
| [`CCD_Oct25_system260522_publication.cys`](CCD_Oct25_system260522_publication.cys) | The Cytoscape session — every systems map, HiDef result, and cosine-similarity subcommunity layout for all 3 cross-species and 4 single-species networks, in one file |

## Inputs (from NDEx, `.cx` format)

The `.cys` session above was built by importing these `.cx` files, deposited at
[`data/05_network/cytoscape/`](../../../data/05_network/cytoscape) — the reproducible inputs to
this manual step. The HiDeF community network and cosine-similarity subgraph they contain come
from real code in `1_NetColoc`'s own source notebook, not something computed within Cytoscape
itself — see [`1_NetColoc/README.md`](../1_NetColoc/README.md)'s Notes for the exact cells and why
it isn't re-executed live here or anywhere else in this repo:

| File pattern | Description |
|---|---|
| `<prefix>_NetColoc_subgraph_CosSim95.cx` | Cosine-similarity subcommunity network |
| `<prefix>_NetColoc_subgraph.cx` | HiDef community network |
| `<prefix>_systems_map.cx` (or `.txt` for `SCH_CCD_Oct25`) | Hierarchical systems map |

| Analysis | Prefix |
|---|---|
| OCD/dogCD | `OCDv4_CCD_Oct25` |
| depression/dogCD | `DEP_CCD_Oct25` |
| schizophrenia/dogCD | `SCH_CCD_Oct25` |
| OCD | `OCD` |
| depression | `DEP` |
| schizophrenia | `SCH` |
| dogCD | `CCD` |

## What happens here

Import each pairing's `.cx` files into Cytoscape, arrange the layout, style nodes/edges/labels by
hand. Final figure polish (fonts, panel assembly, annotations) is done afterward in Adobe
Illustrator (v26.0.2 originally). Neither step is automated or intended to be — the reproducible
artifacts are the `.cx` inputs above, not the hand-laid-out session.
