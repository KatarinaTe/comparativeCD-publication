#!/usr/bin/env python3
"""go_term_comparison.py — GO-enrichment "run1 vs run2" difference volcano plot, generalized by trait.
Source: 05_network/4_Plotting/GO-term_comparison.ipynb — identical logic across the 3 gene-set
pairs (OCD, DEP, SCH), each comparing a single-species-only hierarchy GO enrichment ("run1")
against the corresponding dogCD/CCD cross-species one ("run2").

xlim/ylim: for OCD, the plot limits are computed dynamically from the data; for DEP/SCH, they are
hardcoded instead (`plt.xlim(-20,70)` / `plt.ylim(0,140)`) — a deliberate visual-comparison choice,
not an oversight — reproduced literally via optional --xlim/--ylim rather than made consistent
across all 3 traits.
"""

import argparse

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description="GO-enrichment run1-vs-run2 difference volcano plot for one trait.")
parser.add_argument("run1_go_enrichment", help="Single-species-only hierarchy GO enrichment TSV (index_col=0)")
parser.add_argument("run2_go_enrichment", help="dogCD/CCD cross-species hierarchy GO enrichment TSV (index_col=0)")
parser.add_argument("trait_id", help="Trait label (OCD/DEP/SCH), used only in log output")
parser.add_argument("run1_only_terms_out", help="Output: terms significant in run1 only (single-column, sorted)")
parser.add_argument("full_table_out", help="Output TSV: term, run1, run2, diff, y_max, in_both, run1_only, run2_only, neither")
parser.add_argument("plot_out", help="Output PDF: difference volcano plot")
parser.add_argument("--xlim", type=float, nargs=2, default=None, metavar=("XMIN", "XMAX"), help="Fixed x-axis limits (default: computed symmetrically from the data)")
parser.add_argument("--ylim", type=float, nargs=2, default=None, metavar=("YMIN", "YMAX"), help="Fixed y-axis limits (default: computed from the data)")
args = parser.parse_args()

df1 = pd.read_csv(args.run1_go_enrichment, sep="\t", index_col=0)
df2 = pd.read_csv(args.run2_go_enrichment, sep="\t", index_col=0)
df1["run"] = "run1"
df2["run"] = "run2"
df = pd.concat([df1, df2], ignore_index=True)
term_col = "name"

df["neg_log10_p"] = -np.log10(df["p_value"])
df_sig = df[df["significant"] == True]
max_pivot = df_sig.groupby([term_col, "run"])["neg_log10_p"].max().reset_index()

pivot_table = max_pivot.pivot(index=term_col, columns="run", values="neg_log10_p").fillna(0)
scatter_df = pivot_table

scatter_df["diff"] = scatter_df["run2"] - scatter_df["run1"]
scatter_df["y_max"] = np.maximum(scatter_df["run1"], scatter_df["run2"])

scatter_df["in_both"] = (scatter_df["run1"] > 5) & (scatter_df["run2"] > 5)
scatter_df["run1_only"] = (scatter_df["run1"] > 5) & (scatter_df["run2"] <= 5)
scatter_df["run2_only"] = (scatter_df["run1"] <= 5) & (scatter_df["run2"] > 5)
scatter_df["neither"] = (scatter_df["run1"] <= 5) & (scatter_df["run2"] <= 5)

run1_only_terms = scatter_df[scatter_df["run1_only"]].index.tolist()
pd.Series(run1_only_terms, name="term").sort_values().to_csv(args.run1_only_terms_out, index=False)

x_vals = scatter_df["diff"]
y_vals = scatter_df["y_max"]

full_table = scatter_df.reset_index().rename(columns={term_col: "term"})
full_table = full_table[
    ["term", "run1", "run2", "diff", "y_max", "in_both", "run1_only", "run2_only", "neither"]
].sort_values("diff", ascending=False)
full_table.to_csv(args.full_table_out, sep="\t", index=False)

plt.figure(figsize=(16, 12))

mask_neither = scatter_df["neither"]
plt.scatter(x_vals[mask_neither], y_vals[mask_neither],
            c="lightgray", s=25, alpha=0.5, label=f"Missing both (n={mask_neither.sum():,})")

mask_r1_only = scatter_df["run1_only"]
plt.scatter(x_vals[mask_r1_only], y_vals[mask_r1_only],
            c="blue", s=90, alpha=0.8, edgecolors="darkblue", linewidth=1,
            label=f"Run 1 only (n={mask_r1_only.sum():,})")

mask_r2_only = scatter_df["run2_only"]
plt.scatter(x_vals[mask_r2_only], y_vals[mask_r2_only],
            c="orange", s=90, alpha=0.8, edgecolors="darkorange", linewidth=1,
            label=f"Run 2 only (n={mask_r2_only.sum():,})")

mask_both = scatter_df["in_both"]
plt.scatter(x_vals[mask_both], y_vals[mask_both],
            c="red", s=130, alpha=0.9, edgecolors="black", linewidth=1.8,
            label=f"BOTH runs (n={mask_both.sum():,})", zorder=10)

plt.axvline(x=0, color="black", linestyle="-", lw=2, alpha=0.8, label="No difference")
plt.axhline(y=5, color="gray", linestyle=":", alpha=0.7, lw=2, label="Sig threshold")

plt.xlabel("Run 2 - Run 1 (-log₁₀(p-value))", fontsize=14, fontweight="bold")
plt.ylabel("MAX -log₁₀(p-value) of dominant run", fontsize=14, fontweight="bold")
plt.title("GO Term Significance Difference Between Runs\n"
          f"{len(scatter_df):,} Total Terms | {scatter_df['in_both'].sum():,} Reproducible",
          fontsize=16, fontweight="bold")

plt.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=10, frameon=True, fancybox=True, shadow=True)
plt.grid(True, alpha=0.3)

if args.xlim is not None:
    plt.xlim(args.xlim[0], args.xlim[1])
else:
    max_diff = max(abs(x_vals.max()), abs(x_vals.min()))
    plt.xlim(-max_diff * 1.1, max_diff * 1.1)

if args.ylim is not None:
    plt.ylim(args.ylim[0], args.ylim[1])
else:
    plt.ylim(0, y_vals.max() * 1.05)

plt.tight_layout()
plt.savefig(args.plot_out, format="pdf", bbox_inches="tight", dpi=300)
print(f"Difference volcano saved: {args.plot_out}")

print(f"\nTop 10 by Run 2 - Run 1 difference ({args.trait_id}):")
print(scatter_df.nlargest(10, "diff")[["run1", "run2", "diff", "y_max"]].round(3))
print(f"\nTop 10 by Run 1 - Run 2 difference ({args.trait_id}):")
print(scatter_df.nsmallest(10, "diff")[["run1", "run2", "diff", "y_max"]].round(3))
