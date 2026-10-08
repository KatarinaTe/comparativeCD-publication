#!/usr/bin/env python3
"""control_analysis_plot.py — control-comparisons bar plot (Fig. 1i).
Source: 2.1_OCD_dogCD_Network_Colocalization_260521.ipynb

Consumes the already-deposited `controlanalyses_260521.csv` (data/05_network/) — the
12-combination control loop that produces it calls `calculate_expected_overlap` with an unseeded
permutation (num_reps=1000, no overlap_control), so it is not reproducible; this script reloads
the frozen CSV rather than re-deriving it.

The deposited controlanalyses_260521.csv is semicolon-delimited, read here with `sep=";"`.

Layout matches the published Fig. 1i (2026-10-06): horizontal bars, top to bottom, with
significance stars per the figure legend (* p<2.3x10-5, ** p<3.2x10-14, *** p<1.4x10-68). The
legend thresholds are the rounded p-values of the bars they mark, so p is compared after rounding
to 2 significant figures.
"""

import argparse

import matplotlib
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description="Control-comparisons bar plot (Fig. 1i).")
parser.add_argument("control_results", help="controlanalyses_260521.csv (data/05_network/)")
parser.add_argument("plot_out", help="Output PDF")
args = parser.parse_args()

control_results = pd.read_csv(args.control_results, sep=";", index_col=0)

# (row in CSV, label, colour), top to bottom
rows = [
    ("dCCD-hOCD", "dogCD & OCD", "#ff7100"),
    ("dHeight-hHeight2", "dogSize & height", "#65b400"),
    ("dCCD-hSCH", "dogCD & schizophrenia", "#006a0e"),
    ("hOCD-hSCH", "OCD & schizophrenia", "#006a0e"),
    ("dCCD-hDEP", "dogCD & depression", "#006a0e"),
    ("hOCD-hDEP", "OCD & depression", "#006a0e"),
    ("dCCD-hRA", "dogCD & rheum. arthritis", "#aaaaaa"),
    ("hOCD-hRA", "OCD & rheum. arthritis", "#aaaaaa"),
    ("dCCD-hHeight2", "dogCD & height", "#aaaaaa"),
    ("hOCD-hHeight2", "OCD & height", "#bbbbbb"),
    ("dCCD-hDIA", "dogCD & diabetes", "#bbbbbb"),
    ("hOCD-hDIA", "OCD & diabetes", "#bbbbbb"),
]
main_fig = control_results.loc[[r[0] for r in rows]]
labels = [r[1] for r in rows]
colours = [r[2] for r in rows]


def stars(p):
    p = float(f"{p:.1e}")
    if p <= 1.4e-68:
        return "***"
    if p <= 3.2e-14:
        return "**"
    if p <= 2.3e-5:
        return "*"
    return ""


plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
                     "font.size": 6, "pdf.fonttype": 42, "axes.linewidth": 0.5})
fig, ax = plt.subplots(figsize=(2.6, 2.2))
y = range(len(rows))
ax.barh(
    y=y, left=1, width=main_fig.Mean - 1, height=0.7, color=colours, edgecolor="black", linewidth=0.4,
    xerr=[main_fig.Mean - main_fig.Lower, main_fig.Upper - main_fig.Mean],
    error_kw={"elinewidth": 0.7, "capsize": 2, "capthick": 0.5},
)
for yi, (upper, p) in enumerate(zip(main_fig.Upper, main_fig.p)):
    ax.text(upper + 0.05, yi, stars(p), va="center", ha="left", fontsize=7)
ax.axvline(1, linestyle=(0, (1.5, 1.5)), color="black", linewidth=0.5)
ax.set_yticks(list(y), labels)
ax.tick_params(axis="y", length=0)
ax.invert_yaxis()
ax.set_xlim(0.75, 3.0)
ax.set_xticks([1.0, 1.5, 2.0, 2.5, 3.0])
ax.set_xlabel("Observed/expected size")
for side in ("top", "right", "left"):
    ax.spines[side].set_visible(False)
fig.tight_layout()

fig.savefig(args.plot_out)
