#!/usr/bin/env python3
"""build_heats_matrix.py — builds the normalized adjacency matrix (w_prime) and the individual
heats matrix (w_double_prime) from a fetched interactome.
Source: 1.1_OCD_dogCD_NetColoc_analysis_260521.ipynb

Recomputed live here rather than requiring the deposited 2.8GB matrices: both
get_normalized_adjacency_matrix/get_individual_heats_matrix are pure deterministic linear algebra
(no randomness anywhere). The only cost is a dense O(N^3) matrix inversion (~7x10^12 FLOPs for
N=19,267), a few minutes on one HPC node with a decent BLAS backend, ~10-15GB peak RAM.
Recomputing avoids depending on an external 2.8GB deposit at all.
"""

import argparse
import pickle

import numpy as np
from netcoloc.netprop import get_individual_heats_matrix, get_normalized_adjacency_matrix

parser = argparse.ArgumentParser(description="Build w_prime and the individual heats matrix from a fetched interactome.")
parser.add_argument("interactome_file", help="Pickled networkx.Graph, from fetch_pcnet2.py")
parser.add_argument("--alpha", type=float, default=0.5, help="Heat dissipation coefficient (default: 0.5, matching the notebook's own call)")
parser.add_argument("w_prime_out", help="Output .npy file for the normalized adjacency matrix")
parser.add_argument("indiv_heats_out", help="Output .npy file for the individual heats matrix")
args = parser.parse_args()

with open(args.interactome_file, "rb") as f:
    interactome = pickle.load(f)

print("Calculating w_prime")
w_prime = get_normalized_adjacency_matrix(interactome, conserve_heat=True)

print("Calculating individual heats matrix")
individual_heats_matrix = get_individual_heats_matrix(w_prime, args.alpha)

np.save(args.w_prime_out, w_prime)
np.save(args.indiv_heats_out, individual_heats_matrix)
