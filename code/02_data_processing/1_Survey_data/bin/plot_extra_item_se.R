#!/usr/bin/env Rscript
# plot_extra_item_se.R — F1's extra diagnostic scatter plot (an item's response vs. factor SE)
# Source: 02_data_processing/1_Survey_data/Darwin_dog-CD_3-factorsolution.R, lines 347-348
#
# Only F1 has this plot in the original — invoked for F1 alone (see the main workflow), reading
# the factor's already-saved CCD3F.txt rather than recomputing anything. The preceding
# `table(F1_CCD3F$F1_SE, exclude=NULL)` diagnostic print is not reproduced (see the pipeline
# README's per-factor asymmetries note) — the plot itself is the actual result.

args <- commandArgs(trailingOnly = TRUE)
scores_path <- args[1]
factor_id   <- args[2]
item_col    <- args[3]

out <- read.table(scores_path, header = TRUE, sep = "\t")

png(paste0(factor_id, "_extra_item_se_plot.png"), width = 7, height = 7, units = "in", res = 150)
plot(out[[paste0(factor_id, "_SE")]], out[[paste0("item", item_col)]],
     ylab = paste0("item", item_col), pch = 19, xlab = "SE")
dev.off()
