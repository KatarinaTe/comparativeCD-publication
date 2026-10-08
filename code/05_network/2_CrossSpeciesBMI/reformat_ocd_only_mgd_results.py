#!/usr/bin/env python3
"""reformat_ocd_only_mgd_results.py — normalize the collaborator's re-deposited OCD_ONLY MGD result.
The original data/05_network/OCD_ONLY_hierarchy_full_MGD_enrichment_results.tsv had unrelated LDSC
annotation content (flagged separately). The collaborator re-deposited the real data as
fixedOCD_ONLY_hierarchy_full_MGD_enrichment_results.csv (main branch commit 34ccf40), but in a
different export format than its 5 sibling MGD files: UTF-8 with a BOM, semicolon-delimited,
comma-decimal numbers (e.g. "1749,45"), CRLF line endings, and "TRUE"/"FALSE" instead of
Python-style "True"/"False". Confirmed correct in substance first, not just reformatted blindly:
60 unique communities, matching OCD_ONLY_systems_map.txt's own community count exactly (60 data
rows) — the previous broken file had none of this correspondence.

Run once to produce the properly-formatted TSV (tab-delimited, period-decimal, matching every
other *_hierarchy_full_MGD_enrichment_results.tsv file's exact style) at the same path/filename the
pipeline and the original notebooks expect, replacing the previously-broken content there.
"""

import argparse

import pandas as pd

parser = argparse.ArgumentParser(description="Normalize the re-deposited OCD_ONLY MGD result to match its siblings' format.")
parser.add_argument("input_csv", help="fixedOCD_ONLY_hierarchy_full_MGD_enrichment_results.csv (semicolon-delimited, comma-decimal, UTF-8 BOM)")
parser.add_argument("output_tsv", help="Output TSV, matching the sibling *_hierarchy_full_MGD_enrichment_results.tsv files' format exactly")
args = parser.parse_args()

df = pd.read_csv(args.input_csv, sep=";", index_col=0, encoding="utf-8-sig")
# pandas already parses "TRUE"/"FALSE" into native bool for sig_5e6 during CSV parsing — no
# conversion needed there, only the comma-decimal numeric columns need fixing.

numeric_cols = ["OR", "OR_p", "OR_CI_lower", "OR_CI_upper", "hyper_p"]
for col in numeric_cols:
    df[col] = df[col].astype(str).str.replace(",", ".", regex=False).astype(float)

print(f"{len(df)} rows, {df['name'].nunique()} unique communities")
df.to_csv(args.output_tsv, sep="\t")
