#!/usr/bin/env python3
"""select_igv_individuals.py — for each lead SNP, picks up to N individuals per genotype class
(0/1/2 copies of the counted allele) for manual IGV read-pileup inspection.

Picking multiple individuals per class (not just one) is deliberate: not every dog with an array
genotype call at a lead SNP necessarily has usable BAM read depth AT THAT EXACT POSITION —
sequencing coverage varies — so 3 candidates per class gives a much better chance that at least one
actually shows readable pileup in IGV, without requiring the human to manually scan the full ~1200-
dog cohort themselves.

Input: a plink `--recode A` output (.raw). `--recode A` gives one column per SNP, valued 0/1/2
(copies of the "counted" allele — plink's minor allele by default) or NA (no confident call).
Genotype meaning (which allele is 0 vs 2) isn't resolved here — the human doing the IGV check reads
that directly off the actual base calls in the pileup, cross-referenced against the sumstats A1/A2
columns if needed.

For a SNP with zero individuals in a genotype class, one placeholder row (empty IID) is still
written, so a reader of the output table sees "this class had nobody available" rather than a
silently missing row that could be mistaken for such a check never having run.

Re-running for a locus where none of the first batch showed usable reads: pass --round 2 (3, 4,
...). Each genotype class's full candidate pool is shuffled once, deterministically, from --seed —
round 1 takes the first --n-per-class of that fixed shuffle, round 2 takes the next
--n-per-class-sized slice, and so on, so later rounds never repeat an already-tried individual
without the investigator having to track or pass in an exclude list by hand. Every round is still
fully reproducible from the same (--seed, --round) pair.
"""

import argparse
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser(description="Select up to N individuals per genotype class per lead SNP for IGV review.")
parser.add_argument("raw_file", help="plink --recode A output (.raw)")
parser.add_argument("out_file", help="Output selection table (tsv)")
parser.add_argument("--n-per-class", type=int, default=3, help="Individuals to select per genotype class (default: 3)")
parser.add_argument("--seed", type=int, default=42, help="Random seed for reproducible sampling (default: 42)")
parser.add_argument("--round", type=int, default=1, help="Selection round (1-indexed, default: 1) — increment to get a fresh, non-overlapping batch if an earlier round's picks had no usable reads at this position")
args = parser.parse_args()

if args.round < 1:
    parser.error("--round must be >= 1")

df = pd.read_csv(args.raw_file, sep=r"\s+")
meta_cols = {"FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"}
snp_cols = [c for c in df.columns if c not in meta_cols]

rng = np.random.default_rng(args.seed)

start = (args.round - 1) * args.n_per_class
end = start + args.n_per_class

rows = []
for snp_col in snp_cols:
    for genotype in (0, 1, 2):
        candidates = df.loc[df[snp_col] == genotype, ["FID", "IID"]]
        n_available = len(candidates)
        # Shuffle the FULL candidate pool once (fixed by --seed and loop order, identical across
        # rounds), then slice out this round's window — guarantees rounds never overlap.
        shuffled = candidates.sample(frac=1, random_state=rng)
        batch = shuffled.iloc[start:end]
        if len(batch) > 0:
            for _, row in batch.iterrows():
                rows.append({"SNP": snp_col, "genotype": genotype, "FID": row["FID"], "IID": row["IID"], "n_available_in_class": n_available})
        else:
            rows.append({"SNP": snp_col, "genotype": genotype, "FID": pd.NA, "IID": pd.NA, "n_available_in_class": n_available})

out_df = pd.DataFrame(rows, columns=["SNP", "genotype", "FID", "IID", "n_available_in_class"])
out_df.to_csv(args.out_file, sep="\t", index=False)

print(f"Round {args.round}: selected up to {args.n_per_class} individual(s) per genotype class for {len(snp_cols)} SNP(s).")
for snp_col in snp_cols:
    counts = df[snp_col].value_counts(dropna=True).sort_index().to_dict()
    print(f"  {snp_col}: available per class (0/1/2) = {counts}")
