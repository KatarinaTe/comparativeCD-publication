#!/usr/bin/env python3
"""mgd_significant_communities.py — derive significantly-MGD-enriched communities from a frozen MGD result.
Source: 2.2_OCD_dogCD_Systems_Map_260521.ipynb / 2.3_OCD_ONLY_Systems_Map_260525.ipynb

Reloads the deposited MGD enrichment result rather than re-running `community_term_enrichment`,
which additionally needs the full MGI/MPO/ontology machinery this pipeline doesn't otherwise
build.

Skips the `.merge(final_annotations, ...)` step present in both notebooks: `final_annotations` is
never defined anywhere in either notebook. The "description" column both notebooks report via
`value_counts("description")` is already present in the raw MGD enrichment file itself (one row
per (community, MP-term) pair, each carrying that MP term's own description, e.g. "abnormal
behavior"), and nothing downstream of the merge touches any column that could only come from
`final_annotations`, so this script uses the base MGD file's own "description" column throughout
and omits the merge.

2 of the 6 deposited MGD files (SCH_CCD, DEP_CCD) use a comma as the decimal separator in their
numeric columns (e.g. "2,41E-116"), unlike the other 4 (period-decimal). Read naively, this
silently loads OR_p/OR_CI_lower etc. as strings, which both crashes fdrcorrection() AND makes
sort_values(by="OR_CI_lower") do a lexicographic string sort instead of a numeric one. Coerced
explicitly below.
"""

import argparse

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import fdrcorrection

parser = argparse.ArgumentParser(description="Derive significantly-MGD-enriched communities from a frozen MGD enrichment result.")
parser.add_argument("mgd_enrichment", help="Frozen MGD enrichment TSV (one row per community x MP-term, columns include OR_p, OR_CI_lower, description, name)")
parser.add_argument("comm_results_sign_out", help="Output TSV: sorted by OR_CI_lower, with q (FDR-corrected p), -log10p, and annotate columns added")
args = parser.parse_args()

bmi_pheno_results = pd.read_csv(args.mgd_enrichment, sep="\t", index_col=0)

# Coerce comma-decimal numeric columns (SCH_CCD/DEP_CCD) to floats; a no-op for files that are
# already proper floats.
numeric_cols = ["observed", "total", "OR", "OR_p", "OR_CI_lower", "OR_CI_upper", "hyper_p"]
for col in numeric_cols:
    if bmi_pheno_results[col].dtype == object:
        bmi_pheno_results[col] = bmi_pheno_results[col].str.replace(",", ".", regex=False).astype(float)

print("Description value counts (all MP terms tested):")
print(bmi_pheno_results.value_counts("description"))

comm_results_sign = bmi_pheno_results.sort_values(by="OR_CI_lower", ascending=False)
comm_results_sign["q"] = fdrcorrection(list(comm_results_sign.OR_p.values), method="poscorr")[1]
comm_results_sign["-log10p"] = comm_results_sign.OR_p.apply(lambda x: -1 * np.log10(x))
comm_results_sign["annotate"] = comm_results_sign.apply(lambda x: x["name"] if (x["OR_p"] < 0.05) else "Other", axis=1)

comm_results_sign.to_csv(args.comm_results_sign_out, sep="\t")

print("\nDescription value counts (sorted/annotated result):")
print(comm_results_sign.value_counts("description"))
