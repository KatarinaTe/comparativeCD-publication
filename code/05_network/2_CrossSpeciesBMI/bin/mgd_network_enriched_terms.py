#!/usr/bin/env python3
"""mgd_network_enriched_terms.py -- MP terms enriched anywhere near the top of a systems-map
hierarchy, with per-top-level-community traceability.

Built 2026-10-04 to check the manuscript Results sentence: "The vast majority of MP terms
enriched in the OCD/dogCD combined network (139 of 142) are associated with C185 ... Only 10
terms were for non-neural phenotypes". Neither the frozen MGD enrichment file nor its source
notebook (2.2_OCD_dogCD_Systems_Map_260521.ipynb) contains code that reproduces that exact
139/142/10 breakdown -- the notebook only loads the per-community enrichment table and adds a
`q`/`annotate` column (see mgd_significant_communities.py), nothing that filters to a network-wide
total or classifies terms by community of origin. That tally was assembled manually when the
figure was built (see code/05_network/2_CrossSpeciesBMI/README.md). This script is the closest
reproducible approximation, worked out interactively against the deposited OCD/dogCD file:

  "Enriched in the network" = significant (sig_5e6, i.e. hyper_p < 5e-6) in the union of the
  hierarchy's root community and its immediate child communities (root's direct, non-nested
  sub-communities -- found here as the maximal proper gene-subsets of the root; a community
  nested inside another candidate is excluded as not "immediate"). "Associated with <community>"
  = also significant in that specific community.

On the real OCD/dogCD deposited data this gives 141 terms in the union (root C184 + children
C185/C186/C197), of which 133 are also significant in C185 -- close to but not exactly the
manuscript's 142/139 (the MGD/MPO reference database isn't version-pinned anywhere in this repo,
so a value computed today isn't guaranteed to match one computed in 2025 against the live
database at that time; see data/05_network/README.md). The term *content* matches far more
precisely: of the 122 real (non-degenerate) terms significant in C185, exactly 9 are non-neural
growth/body-size/lethality phenotypes -- versus the manuscript's "only 10" -- and all six named
domains (synaptic physiology, brain morphology, locomotor behaviour, learning/memory, anxiety,
social behaviour) are directly present among C185's top hits.

"Degenerate" terms (flagged, not dropped): MGI/MPO terms annotated to fewer than 10 genes overall
produce unstable, often exactly-zero hyper_p values from a handful of counts (e.g. a term with 1
total annotated gene, observed in this community, trivially is "significant"). Flagged via
`total < 10` so they can be filtered in the output table rather than silently dropped --
consistent with this analysis's own stated background filter (10-2,000 genes per term, see the
manuscript Methods).

Usage: mgd_network_enriched_terms.py <mgd_enrichment.tsv> <hierarchy.tsv> <out.tsv>
  mgd_enrichment.tsv: one row per (community x MP-term), columns include name/description/
                      hyper_p/sig_5e6/total (e.g. 251023_KT_hierarchy_full_MGD_enrichment_results.tsv)
  hierarchy.tsv: community membership, columns include name/CD_MemberList
                 (e.g. CompulsiveNetwork_hierarchy_data_251022.tsv)
"""

import argparse

import pandas as pd

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("mgd_enrichment", help="Per-community MGD/MPO enrichment TSV")
parser.add_argument("hierarchy", help="Hierarchy community-membership TSV (name, CD_MemberList)")
parser.add_argument("out_tsv", help="Output TSV: one row per MP term enriched anywhere in the top-level union")
args = parser.parse_args()

# ---- Find the root community and its immediate children ------------------------------------
hier = pd.read_csv(args.hierarchy, sep="\t", index_col=0)
members = {row["name"]: set(row["CD_MemberList"].split(" ")) for _, row in hier.iterrows()}
sizes = {name: len(genes) for name, genes in members.items()}

root = max(sizes, key=sizes.get)
candidates = {c for c in members if c != root and members[c] < members[root]}
immediate_children = sorted(
    (c for c in candidates if not any(members[c] < members[other] for other in candidates if other != c)),
    key=lambda c: -sizes[c],
)
top_level = [root] + immediate_children
print(f"Root community: {root} ({sizes[root]} genes)")
print(f"Immediate children: {', '.join(f'{c} ({sizes[c]} genes)' for c in immediate_children)}")

# ---- Load enrichment, coercing comma-decimal columns (same issue as mgd_significant_communities.py) --
df = pd.read_csv(args.mgd_enrichment, sep="\t", index_col=0)
df.index.name = "MP_id"
df = df.reset_index()
numeric_cols = ["observed", "total", "OR", "OR_p", "OR_CI_lower", "OR_CI_upper", "hyper_p"]
for col in numeric_cols:
    if df[col].dtype == object:
        df[col] = df[col].str.replace(",", ".", regex=False).astype(float)

# ---- Per-top-level-community columns, one MP term per row -----------------------------------
# Every original column is carried over per community (prefixed by community name), not just
# hyper_p/sig_5e6 -- observed/OR/OR_p/OR_CI_lower/OR_CI_upper/size included too, so nothing from
# the source file is lost and the full per-community stats are directly inspectable alongside the
# trace columns.
original_cols = ["observed", "total", "OR", "OR_p", "OR_CI_lower", "OR_CI_upper", "hyper_p", "sig_5e6", "size"]

descriptions = df.drop_duplicates("MP_id").set_index("MP_id")["description"]
terms = pd.DataFrame(index=sorted(descriptions.index))
terms["description"] = descriptions

sig_cols = []
for community in top_level:
    sub = df[df.name == community].set_index("MP_id")
    for col in original_cols:
        terms[f"{community}_{col}"] = sub[col].reindex(terms.index)
    sig_col = f"{community}_sig_5e6"
    terms[sig_col] = (terms[sig_col] == True)  # noqa: E712 -- also turns NaN (not tested) into False
    sig_cols.append(sig_col)

# A term's background size (`total`) is a property of the MP term itself, constant across
# communities -- root's own `<root>_total` column already carries it (every root-community row is
# present by construction, since top_level's first entry is always the root). Used only to flag
# degenerate terms below, not as a separate new column (that would just duplicate `<root>_total`).
root_total = terms[f"{root}_total"]
terms["degenerate_small_term"] = root_total < 10

terms["n_top_level_sig"] = terms[sig_cols].sum(axis=1)
enriched_in_network = terms[terms["n_top_level_sig"] > 0].copy()
enriched_in_network = enriched_in_network.sort_values(
    by=[f"{c}_hyper_p" for c in top_level]
)

print(f"\n'Enriched in the network' (sig_5e6 in {root} or any immediate child): {len(enriched_in_network)} terms")
for community in top_level:
    n = enriched_in_network[f"{community}_sig_5e6"].sum()
    print(f"  ...of which significant in {community}: {n}")
n_degenerate = enriched_in_network["degenerate_small_term"].sum()
print(f"  ...flagged degenerate ({root}_total < 10): {n_degenerate}")

enriched_in_network.to_csv(args.out_tsv, sep="\t", index_label="MP_id")
print(f"\nWrote {args.out_tsv}")
