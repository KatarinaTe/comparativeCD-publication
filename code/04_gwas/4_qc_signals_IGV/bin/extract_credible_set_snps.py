#!/usr/bin/env python3
"""extract_credible_set_snps.py — pulls the lead-SNP list for IGV review directly out of
3_finemap_susie's own FINEMAP_REGION outputs, rather than waiting on a separate region-selection
step.

PolyFun's finemapper.py (SuSiE method) already tags every SNP with a `CREDIBLE_SET` column in its
own output (0 = not in any credible set; 1, 2, ... = which credible set it belongs to), built
directly from the fitted susieR object's own `sets` attribute. So "which SNPs are worth checking
in IGV" doesn't need the still-on-hold 5_gwas_regions collaborator scripts at all for this
narrower purpose — it's already sitting in every one of 3_finemap_susie's ~104 result files.

Takes any number of finemap_results/*.gz files (3_finemap_susie's own `finemap_results` published
output). One SNP is enough per credible set for an IGV signal check — not every SNP in it — so
within each (file, CREDIBLE_SET) group this keeps only the lowest-P (original GWAS p-value) row,
i.e. the credible set's top/lead SNP. `P` survives unrenamed into the finemap output:
`munge_polyfun_sumstats.py` normalizes whatever p-value column name it finds to exactly `P`, and
finemapper.py's SuSiE branch copies the full sumstats row through into its own output, so `P` is
just carried along. Writes the deduplicated set of lead-SNP IDs (already in chr:pos form, matching
this project's variant-ID convention throughout) — one per line, in plink `--extract` format.
"""

import argparse
import pandas as pd

parser = argparse.ArgumentParser(description="Extract each credible set's top-P SNP from finemapper.py result files.")
parser.add_argument("finemap_results", nargs="+", help="One or more finemapper.py output files (.gz, tab-separated)")
parser.add_argument("out_file", help="Output SNP list, one ID per line (plink --extract format)")
args = parser.parse_args()

snp_ids = set()
for result_file in args.finemap_results:
    df = pd.read_csv(result_file, sep="\t")
    in_credible_set = df.loc[df["CREDIBLE_SET"] > 0]
    top_snps = in_credible_set.loc[in_credible_set.groupby("CREDIBLE_SET")["P"].idxmin(), "SNP"]
    snp_ids.update(top_snps.tolist())

with open(args.out_file, "w") as f:
    for snp_id in sorted(snp_ids):
        f.write(f"{snp_id}\n")

print(f"Found {len(snp_ids)} unique credible-set lead SNP(s) (top P per set) across {len(args.finemap_results)} finemap result file(s).")
