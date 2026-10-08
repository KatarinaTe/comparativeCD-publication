#!/usr/bin/env python3
"""compute_zscores.py — per-trait network propagation z-scores.
Source: 1.1_OCD_dogCD_NetColoc_analysis_260521.ipynb — genericized across all 7 traits by
argument rather than handled per pair.

Seed genes are loaded via pandas (a single-column CSV with a header named "gene"), NOT the plain
whitespace-split text file netcoloc's own top-level netprop_zscore() wrapper expects — this calls
calculate_heat_zscores() directly and does its own seed-gene loading and interactome intersection
first (redundant with calculate_heat_zscores()'s own internal np.intersect1d call, but kept for
fidelity).

alpha=0.5 and random_seed=1 are netcoloc's own function defaults, passed explicitly here as script
arguments (still defaulting to the same values) so the exact values used don't silently depend on
whatever netcoloc's own defaults happen to be in a future release.

Output matches the deposited historical z-score CSVs: a single 'z' column, written with the
default pandas index (gene names).
"""

import argparse
import pickle

import numpy as np
import pandas as pd
from netcoloc.netprop_zscore import calculate_heat_zscores

parser = argparse.ArgumentParser(description="Compute network propagation z-scores for one trait's seed genes.")
parser.add_argument("interactome_file", help="Pickled networkx.Graph, from fetch_pcnet2.py")
parser.add_argument("indiv_heats_matrix_file", help=".npy individual heats matrix, from build_heats_matrix.py")
parser.add_argument("seed_gene_file", help="Single-column CSV with a 'gene' header (e.g. OCDgenes_v4.txt)")
parser.add_argument("--num-reps", type=int, default=1000, help="Null model repetitions (default: 1000, matching the notebook)")
parser.add_argument("--minimum-bin-size", type=int, default=100, help="Degree-matching bin size (default: 100, matching the notebook)")
parser.add_argument("--alpha", type=float, default=0.5, help="Heat dissipation coefficient (default: 0.5, netcoloc's own default)")
parser.add_argument("--random-seed", type=int, default=1, help="RNG seed (default: 1, netcoloc's own default, never overridden in the notebook)")
parser.add_argument("z_scores_out", help="Output CSV (single 'z' column, gene-name index)")
args = parser.parse_args()

with open(args.interactome_file, "rb") as f:
    interactome = pickle.load(f)
int_nodes = list(interactome.nodes)
degrees = dict(interactome.degree)

individual_heats_matrix = np.load(args.indiv_heats_matrix_file)

seed_df = pd.read_csv(args.seed_gene_file)
seed_df.index = seed_df["gene"]
seed_genes = seed_df.index.tolist()
print(f"Number of seed genes: {len(seed_genes)}")
seed_genes = list(np.intersect1d(seed_genes, int_nodes))
print(f"Number of seed genes in interactome: {len(seed_genes)}")

print("Calculating z-scores")
z_scores, final_heat, random_final_heats = calculate_heat_zscores(
    individual_heats_matrix,
    int_nodes,
    degrees,
    seed_genes,
    num_reps=args.num_reps,
    alpha=args.alpha,
    minimum_bin_size=args.minimum_bin_size,
    random_seed=args.random_seed,
)

z_scores_df = pd.DataFrame({"z": z_scores})
z_scores_df.to_csv(args.z_scores_out)
