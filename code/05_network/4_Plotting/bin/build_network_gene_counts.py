#!/usr/bin/env python3
"""build_network_gene_counts.py -- live-computed replacement for Fig. 2a's hardcoded
captured/total gene-count literals in plot_network_overlap.R.

Source: this repo's own verification work (REPRODUCIBILITY_AUDIT.md, 2026-09-27 entries) found
that one of the 7 hardcoded (captured, total) pairs was wrong -- Depression/dogCD's "total" was
coded as 387 but a naive z_comb>3 recount from the deposited data gave 495. Tracing this down (via
the real source notebooks pulled from Pelle) found the true cause: Depression/dogCD's network was
deliberately built at z_comb=4 ("zthresh=4 # default = 3" in DEP_CCD_Oct25_NetColoc_analysis_251022.ipynb),
not the z_comb=3 used everywhere else -- 387 was correct, but only at the right threshold. A second,
independent bug was also found: the per-network seed-gene breakdown (dog seeds in network) had been
copied as a constant "8" across all three cross-species rows in a hand-maintained tracking sheet,
when the real values are 8/23/16 for OCD/Depression/Schizophrenia respectively.

This script recomputes all 7 networks directly from the same deposited z-score, seed-gene, and
systems-map files used throughout 05_network, so a future change to any of those inputs is
reflected automatically instead of needing another hand-verification pass.
"""

import argparse

import pandas as pd

parser = argparse.ArgumentParser(description="Recompute Fig. 2a's per-network gene/seed counts live.")

# Single-species: z-score file, seed-gene file, systems-map (hierarchy) gene list
parser.add_argument("--dogcd-z", required=True)
parser.add_argument("--dogcd-seeds", required=True)
parser.add_argument("--dogcd-seeds-header", action="store_true", help="Seed file has a header row (e.g. 'gene') to skip.")
parser.add_argument("--dogcd-hier", required=True)

parser.add_argument("--ocd-z", required=True)
parser.add_argument("--ocd-seeds", required=True)
parser.add_argument("--ocd-seeds-header", action="store_true")
parser.add_argument("--ocd-hier", required=True)

parser.add_argument("--dep-z", required=True)
parser.add_argument("--dep-seeds", required=True)
parser.add_argument("--dep-seeds-header", action="store_true")
parser.add_argument("--dep-hier", required=True)

parser.add_argument("--sch-z", required=True)
parser.add_argument("--sch-seeds", required=True)
parser.add_argument("--sch-seeds-header", action="store_true")
parser.add_argument("--sch-hier", required=True)

# Cross-species: combined zcomb_z12 file, systems-map (hierarchy) gene list, chosen z_comb threshold
# (Depression/dogCD's default of 4 reflects the real, deliberate choice found in its source
# notebook -- not the same as the other two, which use netcoloc's own default of 3.)
parser.add_argument("--ocd-ccd-zcomb", required=True)
parser.add_argument("--ocd-ccd-hier", required=True)
parser.add_argument("--ocd-ccd-z-comb-threshold", type=float, default=3.0)

parser.add_argument("--dep-ccd-zcomb", required=True)
parser.add_argument("--dep-ccd-hier", required=True)
parser.add_argument("--dep-ccd-z-comb-threshold", type=float, default=4.0)

parser.add_argument("--sch-ccd-zcomb", required=True)
parser.add_argument("--sch-ccd-hier", required=True)
parser.add_argument("--sch-ccd-z-comb-threshold", type=float, default=3.0)

parser.add_argument("--z1-threshold", type=float, default=1.5)
parser.add_argument("--z2-threshold", type=float, default=1.5)
parser.add_argument("--single-species-threshold", type=float, default=3.0)

parser.add_argument("output", help="Output CSV: network_gene_counts.csv")
args = parser.parse_args()


def load_gene_set(path, has_header):
    with open(path) as f:
        lines = [l.strip() for l in f if l.strip()]
    return set(lines[1:] if has_header else lines)


def single_species_row(name, species, z_file, seed_file, seed_has_header, hier_file, threshold):
    z = pd.read_csv(z_file, index_col=0)
    z.columns = ["z"]
    seeds = load_gene_set(seed_file, seed_has_header)
    hier = load_gene_set(hier_file, False)

    n_tot = int((z["z"] > threshold).sum())
    n_hier = len(hier)

    seeds_present = seeds & set(z.index)
    seeds_in_net = {g for g in seeds_present if z.loc[g, "z"] > threshold}
    seeds_in_hier = seeds_in_net & hier

    dog_tot, hum_tot = (len(seeds_present), 0) if species == "dog" else (0, len(seeds_present))
    dog_net, hum_net = (len(seeds_in_net), 0) if species == "dog" else (0, len(seeds_in_net))
    dog_hier, hum_hier = (len(seeds_in_hier), 0) if species == "dog" else (0, len(seeds_in_hier))

    return dict(
        Network=name, Type="single", species=species,
        **{"network thresholds": f"NPS{'d' if species=='dog' else 'h'}>{threshold:g}"},
        **{"genes in network": n_tot, "genes in hierarchy": n_hier},
        **{"seeds total": len(seeds_present), "dog seeds total": dog_tot, "human seeds total": hum_tot},
        **{"seeds in network": len(seeds_in_net), "dog seeds in network": dog_net, "human seeds in network": hum_net},
        **{"dog seeds in hierarchy": dog_hier, "human seeds in hierarchy": hum_hier},
    )


def cross_species_row(name, zcomb_file, hier_file, dog_seeds, human_seeds, threshold, z1_thr, z2_thr):
    d = pd.read_csv(zcomb_file, sep="\t", index_col=0)
    mask = (d["NPS_h"] > z2_thr) & (d["NPS_r"] > z1_thr) & (d["NPS_hr"] > threshold)
    network_genes = set(d.index[mask])
    hier = load_gene_set(hier_file, False)

    dog_present = dog_seeds & set(d.index)
    hum_present = human_seeds & set(d.index)
    dog_in_net = network_genes & dog_present
    hum_in_net = network_genes & hum_present
    dog_in_hier = dog_in_net & hier
    hum_in_hier = hum_in_net & hier

    return dict(
        Network=name, Type="coloc", species="cross",
        **{"network thresholds": f"NPSd>{z1_thr:g}, NPSh>{z2_thr:g}, NPSdh>{threshold:g}"},
        **{"genes in network": len(network_genes), "genes in hierarchy": len(hier)},
        **{"seeds total": len(dog_present) + len(hum_present), "dog seeds total": len(dog_present), "human seeds total": len(hum_present)},
        **{"seeds in network": len(dog_in_net) + len(hum_in_net), "dog seeds in network": len(dog_in_net), "human seeds in network": len(hum_in_net)},
        **{"dog seeds in hierarchy": len(dog_in_hier), "human seeds in hierarchy": len(hum_in_hier)},
    )


dogcd_seeds = load_gene_set(args.dogcd_seeds, args.dogcd_seeds_header)
ocd_seeds = load_gene_set(args.ocd_seeds, args.ocd_seeds_header)
dep_seeds = load_gene_set(args.dep_seeds, args.dep_seeds_header)
sch_seeds = load_gene_set(args.sch_seeds, args.sch_seeds_header)

rows = [
    single_species_row("dogCD", "dog", args.dogcd_z, args.dogcd_seeds, args.dogcd_seeds_header, args.dogcd_hier, args.single_species_threshold),
    single_species_row("OCD", "human", args.ocd_z, args.ocd_seeds, args.ocd_seeds_header, args.ocd_hier, args.single_species_threshold),
    single_species_row("Depression", "human", args.dep_z, args.dep_seeds, args.dep_seeds_header, args.dep_hier, args.single_species_threshold),
    single_species_row("Schizophrenia", "human", args.sch_z, args.sch_seeds, args.sch_seeds_header, args.sch_hier, args.single_species_threshold),
    cross_species_row("OCD/dogCD", args.ocd_ccd_zcomb, args.ocd_ccd_hier, dogcd_seeds, ocd_seeds, args.ocd_ccd_z_comb_threshold, args.z1_threshold, args.z2_threshold),
    cross_species_row("Depression/dogCD", args.dep_ccd_zcomb, args.dep_ccd_hier, dogcd_seeds, dep_seeds, args.dep_ccd_z_comb_threshold, args.z1_threshold, args.z2_threshold),
    cross_species_row("Schizophrenia/dogCD", args.sch_ccd_zcomb, args.sch_ccd_hier, dogcd_seeds, sch_seeds, args.sch_ccd_z_comb_threshold, args.z1_threshold, args.z2_threshold),
]

df = pd.DataFrame(rows)
df["% genes in hierarchy"] = (100 * df["genes in hierarchy"] / df["genes in network"]).round(1).astype(str) + "%"
df["% seed genes in network"] = (100 * df["seeds in network"] / df["seeds total"]).round(1).astype(str) + "%"

df.to_csv(args.output, index=False)
print(df.to_string(index=False))
