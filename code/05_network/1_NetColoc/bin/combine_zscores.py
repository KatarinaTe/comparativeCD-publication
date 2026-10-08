#!/usr/bin/env python3
"""combine_zscores.py — combine dogCD/CCD and one human-trait z-score vector into a zcomb_z12 table.
Source: 1.1_OCD_dogCD_NetColoc_analysis_260521.ipynb — full_join(ccd, ocd, by="gene") then
zcomb = D1_z * D2_z, columns renamed to NPS_r/NPS_h/NPS_hr ("h" = human, "r" = dog, "hr" =
dog-human; "r" is a historical holdover from an earlier version of this pipeline that compared
against rat, kept for column-naming consistency across all 3 deposited *_zcomb_z12* files even
though every non-human trait here is the dog CCD set).

Generalized here by argument to cover all 3 trait pairings (OCD/DEP/SCH vs. CCD) rather than
handled per-pair, since CCD is always the "r" side.

R's full_join keeps every gene present in either input (outer join), with NA where one side is
missing; pandas' outer join + multiplication does the same (NaN propagates through '*' exactly like
R's NA). Row order follows the "r" (dog/CCD) input first, matching full_join(ccd, ocd, ...)'s own
argument order.
"""

import argparse

import pandas as pd

parser = argparse.ArgumentParser(description="Combine dogCD/CCD and one human-trait z-score vector into a zcomb_z12 table.")
parser.add_argument("r_zscores", help="Dog/CCD z-scores CSV (index_col=0, single 'z' column), from compute_zscores.py")
parser.add_argument("h_zscores", help="Human-trait z-scores CSV (index_col=0, single 'z' column), from compute_zscores.py")
parser.add_argument("zcomb_out", help="Output TSV: gene, NPS_r, NPS_h, NPS_hr")
args = parser.parse_args()

r_df = pd.read_csv(args.r_zscores, index_col=0).rename(columns={"z": "NPS_r"})
h_df = pd.read_csv(args.h_zscores, index_col=0).rename(columns={"z": "NPS_h"})

nps_df = r_df.join(h_df, how="outer")
nps_df["NPS_hr"] = nps_df["NPS_r"] * nps_df["NPS_h"]

nps_df.index.name = "gene"
nps_df.reset_index().to_csv(args.zcomb_out, sep="\t", index=False)
