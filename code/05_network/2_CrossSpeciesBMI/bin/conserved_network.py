#!/usr/bin/env python3
"""conserved_network.py — extract the dogCD/CCD-OCD conserved network subgraph, plus Fig. 1e's
threshold-passing gene counts and overlap significance test.
Source: 2.1_OCD_dogCD_Network_Colocalization_260521.ipynb

Reproduces the notebook's literal manual threshold (z_human > 1.5 & z_rat > 1.5 &
z_human*z_rat > 3) directly, rather than calling CrossSpeciesBMI's own
`calculate_network_overlap` helper — mathematically equivalent here since both individual
thresholds already force z_human/z_rat positive, making that helper's extra `*(z1>0)*(z2>0)`
masking a no-op for this specific threshold combination.

Also computes the gene counts and hypergeometric overlap test behind Fig. 1e ("genes passing NPS
thresholds... NPSh > 1.5 (human, blue), NPSd > 1.5 (dog, red), and the colocalized network
(orange)... Overlap p-values were calculated via hypergeometric test") — real code
(`scipy.stats.hypergeom`) the source notebook imports but never actually calls anywhere; no other
script in this repo computes it either. Population size for the test is every gene that received a
network-propagation z-score (`len(nps_df)`), i.e. every PCNet2.0 node reachable by the propagation.
This produces the counts and the significance value; the panel's visual layout itself (like Fig.
1a/1d/3c/5) is still assembled manually from these numbers.
"""

import argparse
import pickle

import networkx as nx
import pandas as pd
from scipy.stats import hypergeom

parser = argparse.ArgumentParser(description="Extract the dogCD/CCD-OCD conserved network subgraph and Fig. 1e's threshold gene-set counts from a zcomb_z12 table.")
parser.add_argument("zcomb_z12", help="gene, NPS_r (dog), NPS_h (human), NPS_hr (combined) TSV — from 1_NetColoc's COMBINE_ZSCORES")
parser.add_argument("interactome", help="Pickled networkx.Graph (PCNet2.0), from 1_NetColoc's FETCH_PCNET2")
parser.add_argument("edgelist_out", help="Output TSV: conserved-network edgelist (nx.to_pandas_edgelist format)")
parser.add_argument("threshold_counts_out", help="Output TSV: Fig. 1e gene counts + hypergeometric overlap p-value")
args = parser.parse_args()

nps_df = pd.read_csv(args.zcomb_z12, sep="\t", index_col=0)

human_threshold = nps_df["NPS_h"] > 1.5
dog_threshold = nps_df["NPS_r"] > 1.5

conserved_network_genes = nps_df[
    human_threshold & dog_threshold & (nps_df["NPS_h"] * nps_df["NPS_r"] > 3)
].index.values
print(f"Number of conserved network genes: {len(conserved_network_genes)}")

with open(args.interactome, "rb") as f:
    interactome = pickle.load(f)

G_conserved = interactome.subgraph(conserved_network_genes)
print(f"Nodes: {len(G_conserved.nodes())}, Edges: {len(G_conserved.edges())}")

nx.to_pandas_edgelist(G_conserved).to_csv(args.edgelist_out, sep="\t")

#### Fig. 1e: threshold-passing gene counts + overlap significance ####
population_size = len(nps_df)
n_human = int(human_threshold.sum())  # NPSh > 1.5 ("human, blue")
n_dog = int(dog_threshold.sum())      # NPSd > 1.5 ("dog, red")
n_overlap = int((human_threshold & dog_threshold).sum())
n_colocalized = len(conserved_network_genes)  # additionally NPShr > 3.0 ("colocalized network, orange")

# One-sided test: is the observed overlap between the two threshold sets bigger than expected by
# chance, drawing n_dog genes from a population of population_size with n_human "successes"?
overlap_p = hypergeom.sf(n_overlap - 1, population_size, n_human, n_dog)

with open(args.threshold_counts_out, "w") as f:
    f.write(f"population_size\t{population_size}\n")
    f.write(f"human_threshold_NPSh_gt_1.5\t{n_human}\n")
    f.write(f"dog_threshold_NPSd_gt_1.5\t{n_dog}\n")
    f.write(f"overlap_both_thresholds\t{n_overlap}\n")
    f.write(f"colocalized_network_NPShr_gt_3\t{n_colocalized}\n")
    f.write(f"overlap_hypergeometric_p\t{overlap_p}\n")
print(f"Threshold overlap: {n_overlap} genes (hypergeometric p={overlap_p:.3g})")
