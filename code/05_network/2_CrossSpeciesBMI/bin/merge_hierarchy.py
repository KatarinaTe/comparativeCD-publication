#!/usr/bin/env python3
"""merge_hierarchy.py — combine hierarchy community membership with per-community seed-gene fractions.
Source: 2.2_OCD_dogCD_Systems_Map_260521.ipynb — `hier_df_genes` (CD_MemberList, represents,
indexed by community name) from the frozen hierarchy data, joined (inner) against
`frac_d1_seeds`/`frac_d2_seeds` from the systems map.
"""

import argparse

import pandas as pd

parser = argparse.ArgumentParser(description="Combine hierarchy community membership with per-community seed-gene fractions.")
parser.add_argument("hierarchy", help="Hierarchy data TSV (e.g. CompulsiveNetwork_hierarchy_data_251022.tsv), columns include name/CD_MemberList/represents")
parser.add_argument("systems_map", help="Systems map TSV, indexed by community name, columns include frac_d1_seeds/frac_d2_seeds")
parser.add_argument("merged_out", help="Output TSV: CD_MemberList, represents, frac_d1_seeds, frac_d2_seeds, indexed by community name")
args = parser.parse_args()

G_BMI_hier = pd.read_csv(args.hierarchy, sep="\t", index_col=0)
hier_df_genes = G_BMI_hier.loc[:, ["name", "CD_MemberList", "represents"]]
hier_df_genes = hier_df_genes.set_index("name")

G_BMI_hier_df = pd.read_csv(args.systems_map, sep="\t", index_col=0)
G_hier_df_genes = G_BMI_hier_df.loc[:, ["frac_d1_seeds", "frac_d2_seeds"]]

merged_df = pd.merge(hier_df_genes, G_hier_df_genes, left_index=True, right_index=True, how="inner")

merged_df["CD_MemberList"] = merged_df["CD_MemberList"].astype(str)
merged_df["represents"] = merged_df["represents"].astype(str)
merged_df["frac_d1_seeds"] = pd.to_numeric(merged_df["frac_d1_seeds"], errors="coerce")
merged_df["frac_d2_seeds"] = pd.to_numeric(merged_df["frac_d2_seeds"], errors="coerce")

print(f"Merged {len(merged_df)} communities")
merged_df.to_csv(args.merged_out, sep="\t")
