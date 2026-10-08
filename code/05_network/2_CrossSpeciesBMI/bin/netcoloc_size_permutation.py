#!/usr/bin/env python3
"""netcoloc_size_permutation.py — dogCD/CCD-OCD colocalized network size vs. a permuted null (Figure 1g).
Source: 2.1_OCD_dogCD_Network_Colocalization_260521.ipynb — calls
`calculate_expected_overlap(z_rat, z_human, z_score_threshold=3, z1_threshold=1.5, z2_threshold=1.5,
num_reps=10000, overlap_control="bin", seed1=dog_seeds, seed2=human_seeds)` then
`plot_permutation_histogram(permuted, observed, ...)`.

`calculate_expected_overlap` is a CrossSpeciesBMI-specific reimplementation, not netcoloc's own
`network_colocalization` module (which this notebook doesn't even import).

The "observed" network overlap size is fully deterministic (computed once from the real z-scores,
no randomness) and is reproduced exactly. The permuted null distribution is NOT reproducible — the
source shuffles with `random.shuffle()` and never seeds it — so this script generates a genuine,
fresh permutation each run. The resulting p-value will therefore vary slightly run-to-run around
the same true value.
"""

import argparse
import random

import matplotlib
import numpy as np
import pandas as pd
from scipy.stats import norm

matplotlib.use("Agg")
import math

import matplotlib.pyplot as plt
import seaborn as sns

parser = argparse.ArgumentParser(description="dogCD/CCD-OCD colocalized network size vs. a permuted null (Figure 1g).")
parser.add_argument("zcomb_z12", help="gene, NPS_r (dog), NPS_h (human), NPS_hr (combined) TSV — from 1_NetColoc's COMBINE_ZSCORES")
parser.add_argument("dog_seed_genes", help="Dog/CCD seed gene CSV (header 'gene'), e.g. CCDgenes_Oct25.txt")
parser.add_argument("human_seed_genes", help="Human/OCD seed gene CSV (header 'gene'), e.g. OCDgenes_v4.txt")
parser.add_argument("histogram_out", help="Output PDF: observed-vs-permuted histogram (Figure 1g)")
parser.add_argument("observed_out", help="Output text file: observed network overlap size and p-value")
parser.add_argument("--z-score-threshold", type=float, default=3.0)
parser.add_argument("--z1-threshold", type=float, default=1.5)
parser.add_argument("--z2-threshold", type=float, default=1.5)
parser.add_argument("--num-reps", type=int, default=10000)
args = parser.parse_args()


def calculate_network_overlap(z_scores_1, z_scores_2, z_score_threshold, z1_threshold, z2_threshold):
    z_combined = z_scores_1 * z_scores_2 * (z_scores_1 > 0) * (z_scores_2 > 0)
    return z_combined[
        (z_combined >= z_score_threshold) & (z_scores_1 > z1_threshold) & (z_scores_2 > z2_threshold)
    ].index.tolist()


def get_p_from_permutation_results(observed, permuted):
    p = norm.sf((observed - np.mean(permuted)) / np.std(permuted))
    try:
        p = round(p, 4 - int(math.floor(math.log10(abs(p)))) - 1)
    except ValueError:
        print("Cannot round result, p=", p)
    return p


nps_df = pd.read_csv(args.zcomb_z12, sep="\t", index_col=0)
z_rat = nps_df["NPS_r"]
z_human = nps_df["NPS_h"]

dog_seeds = pd.read_csv(args.dog_seed_genes)["gene"].tolist()
human_seeds = pd.read_csv(args.human_seed_genes)["gene"].tolist()

seed_overlap = list(set(dog_seeds).intersection(set(human_seeds)))
print(f"Overlap seed genes: {len(seed_overlap)}")

overlap_z_rat = z_rat.loc[seed_overlap]
overlap_z_human = z_human.loc[seed_overlap]
z_rat = z_rat.drop(seed_overlap)
z_human = z_human.drop(seed_overlap)

observed = len(calculate_network_overlap(z_rat, z_human, args.z_score_threshold, args.z1_threshold, args.z2_threshold))
observed += len(calculate_network_overlap(overlap_z_rat, overlap_z_human, args.z_score_threshold, args.z1_threshold, args.z2_threshold))
print(f"Observed network overlap size: {observed}")

z_rat_arr = np.array(z_rat)
overlap_z_rat_arr = np.array(overlap_z_rat)
permuted = np.zeros(args.num_reps)
for i in range(args.num_reps):
    random.shuffle(z_rat_arr)
    perm_size = len(calculate_network_overlap(pd.Series(z_rat_arr, index=z_rat.index), z_human, args.z_score_threshold, args.z1_threshold, args.z2_threshold))
    random.shuffle(overlap_z_rat_arr)
    perm_size += len(calculate_network_overlap(pd.Series(overlap_z_rat_arr, index=overlap_z_rat.index), overlap_z_human, args.z_score_threshold, args.z1_threshold, args.z2_threshold))
    permuted[i] = perm_size

p_value = get_p_from_permutation_results(observed, permuted)
print(f"p-value: {p_value}")

with open(args.observed_out, "w") as f:
    f.write(f"observed\t{observed}\n")
    f.write(f"p_value\t{p_value}\n")
    f.write(f"permuted_mean\t{np.mean(permuted)}\n")
    f.write(f"permuted_std\t{np.std(permuted)}\n")

plt.figure(figsize=(5, 4))
dfig = sns.histplot(permuted, label="Permuted", alpha=0.4, stat="density", bins=25, kde=True, edgecolor="w", color="dimgrey")
plt.rcParams.update({"mathtext.default": "regular"})
plt.xlabel("Size of colocalized network", fontsize=16)
diff = max(observed, max(permuted)) - min(permuted)
plt.arrow(
    x=observed, y=dfig.dataLim.bounds[3] / 2, dx=0, dy=-1 * dfig.dataLim.bounds[3] / 2, label="Observed",
    width=diff / 100, head_width=diff / 15, head_length=dfig.dataLim.bounds[3] / 20, overhang=0.5,
    length_includes_head=True, color="#F5793A", zorder=50,
)
plt.ylabel("Density", fontsize=16)
plt.legend(fontsize=12, loc=(0.6, 0.75))
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.locator_params(axis="y", nbins=6)
plt.title(f"(p={p_value})", fontsize=16)
plt.tight_layout()
plt.savefig(args.histogram_out)
