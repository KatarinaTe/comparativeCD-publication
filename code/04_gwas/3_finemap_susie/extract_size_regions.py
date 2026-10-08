#!/usr/bin/env python3
"""Regenerates assets/finemap_regions_size.csv from data/04_gwas/SIZE_260603_CSnonfunct.xlsx —
SIZE's own already-run SuSiE finemapping results (added to `main` 2026-09-02, after SIZE had
already been extended through 1_mlma-loco's clumping/export but was still excluded from
3_finemap_susie's harmonize/munge/finemap steps for lack of any region-window source).

The original finemapNOFUNCT_multiple_CCDregions*.sh scripts never covered SIZE at all — no
literal source exists for what windows to finemap it at. This spreadsheet is direct evidence
that SIZE WAS finemapped in the real analysis, and its own already-run region filenames
double as that missing window list: they're formatted "SIZE_nonfunct_<chr>[._]<start>[-_]<end>Mb.gz"
(chr-vs-Mb separator is inconsistently "." or "_" across rows — a hand-named file, not a script
output convention), e.g. "SIZE_nonfunct_9.14.2_15.2Mb.gz" -> chr9, 14.2-15.2Mb.

Deliberately does NOT reuse the credible-set boundary columns ("CS boundaries - start/stop") as
the finemapper.py --start/--end window — those are the CREDIBLE SET's own (much narrower) bounds,
not the region actually analysed to find it. Re-running finemapper.py on just the credible set's
own span would be a materially different, incorrect analysis. The literal window comes from the
region filename itself.

The spreadsheet interleaves two logically separate tables in one sheet: a messy per-SNP raw
finemapper.py-output dump (left, columns A-N, one block per region with an intervening blank row)
and a clean one-row-per-region summary (right, columns R-W: Region, Region size, CS nr SNPs, CS
boundaries start/stop, CS size). This script reads only the clean right-hand block — the left
block's per-SNP detail will be regenerated for real by actually re-running FINEMAP_REGION on these
windows, not imported as static data.

N (sample size) is read directly from the left block's own N column (2330, verified identical
across every SIZE row in the file) — a literal value from the source, matching how the other
104 rows' N values are hardcoded literals from their own original scripts, not computed at
pipeline-run time.

Requires pandas + openpyxl (both in the polyfun container's own env — openpyxl added
specifically for this script; see env/04_gwas/3_finemap_susie/polyfun.yml).
"""

import os
import re

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
XLSX_PATH = os.path.join(HERE, "..", "..", "..", "data", "04_gwas", "SIZE_260603_CSnonfunct.xlsx")
OUT_PATH = os.path.join(HERE, "assets", "finemap_regions_size.csv")

REGION_PATTERN = re.compile(r"SIZE_nonfunct_(\d+)[._](\d+(?:\.\d+)?)[-_](\d+(?:\.\d+)?)Mb\.gz")

df = pd.read_excel(XLSX_PATH, header=None)

# Left block: N is in column 4, wherever a real data row (not the "N" sub-header) has a value.
n_values = df.iloc[:, 4].dropna()
n_values = n_values[n_values != "N"].unique()
if len(n_values) != 1:
    raise ValueError(f"Expected exactly one distinct N value across the file, found: {n_values}")
n_size = int(n_values[0])

# Right block: columns 17-22 = Region, Region size, CS nr SNPs, CS boundaries start, stop, CS size.
region_col = df.iloc[:, 17].dropna()
region_names = [v for v in region_col if isinstance(v, str) and v.startswith("SIZE_nonfunct_")]

rows = []
for region_name in region_names:
    m = REGION_PATTERN.match(region_name)
    if not m:
        raise ValueError(f"Region filename didn't match the expected pattern: {region_name}")
    chrom, start_mb, end_mb = m.groups()
    start_bp = round(float(start_mb) * 1_000_000)
    end_bp = round(float(end_mb) * 1_000_000)
    label = f"{chrom}.{start_mb}-{end_mb}Mb"
    rows.append({
        "pheno_id": "SIZE",
        "chr": chrom,
        "start": start_bp,
        "end": end_bp,
        "n": n_size,
        "label": label,
        "out_base": f"SIZE_NO_PRIORS.{label}.gz",
        "source_file": "SIZE_260603_CSnonfunct.xlsx",
    })

out_df = pd.DataFrame(rows, columns=["pheno_id", "chr", "start", "end", "n", "label", "out_base", "source_file"])
out_df.to_csv(OUT_PATH, index=False)

print(f"Wrote {len(rows)} SIZE region(s) to {OUT_PATH}")
