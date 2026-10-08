#!/usr/bin/env python3
"""Regenerates assets/finemap_regions.csv from the original finemapNOFUNCT_multiple_CCDregions*.sh
scripts in this directory.

Not a pipeline process script (nothing under bin/ calls this) — a one-off provenance/regeneration
tool, checked in so the 104-row manifest that 3_finemap_susie.nf reads can be independently
re-derived and verified against the original shell scripts, rather than trusted as a hand-typed
transcription. Run from this directory: `python3 extract_finemap_regions.py`.
"""

import re
import csv
import glob
import os

HERE = os.path.dirname(os.path.abspath(__file__))

files = sorted(glob.glob(os.path.join(HERE, 'finemapNOFUNCT_multiple_CCDregions*.sh')))
rows = []
for fpath in files:
    text = open(fpath).read()
    calls = re.findall(r'python3\s+finemapper\.py\s+(.*?)(?=\npython3\s+finemapper\.py|\Z)', text, re.S)
    for c in calls:
        def grab(flag):
            m = re.search(r'--' + flag + r'\s+(\S+)', c)
            return m.group(1) if m else None
        geno = grab('geno')
        sumstats = grab('sumstats')
        n = grab('n')
        chrom = grab('chr')
        start = grab('start')
        end = grab('end')
        out = grab('out')
        if not (geno and sumstats and out):
            continue
        pheno = re.search(r'/([A-Za-z0-9]+)_sumstats_munged\.parquet', sumstats).group(1)
        out_base = out.split('/')[-1]
        label = re.sub(r'\.gz$', '', out_base)
        label = re.sub(r'^' + re.escape(pheno) + r'_NO_PRIORS\.', '', label)
        rows.append({
            'pheno_id': pheno,
            'chr': chrom,
            'start': start,
            'end': end,
            'n': n,
            'label': label,
            'out_base': out_base,
            'source_file': os.path.basename(fpath),
        })

out_path = os.path.join(HERE, 'assets', 'finemap_regions.csv')
with open(out_path, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=['pheno_id', 'chr', 'start', 'end', 'n', 'label', 'out_base', 'source_file'])
    w.writeheader()
    w.writerows(rows)

print(f"Wrote {len(rows)} rows to {out_path}")
