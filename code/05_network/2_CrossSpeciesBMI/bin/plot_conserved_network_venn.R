#!/usr/bin/env Rscript
# plot_conserved_network_venn.R -- Fig. 1e: genes passing network proximity score (NPS)
# thresholds after network propagation. A plain 2-circle Venn: NPSh > 1.5 only (human, blue),
# NPSd > 1.5 only (dog, red), and both plus NPShd > 3.0 (colocalized network, orange).
#
# Source: bin/conserved_network.py's own fig1e_threshold_counts.tsv output (added 2026-09-21,
# CONSERVED_NETWORK). Confirmed against the real embedded Fig. 1 image
# (NATURE_Main_Figures.docx, 2026-09-22): the panel shows exactly 1369 (NPSh-only) / 184
# (colocalized) / 1753 (NPSd-only) -- matching human_threshold_NPSh_gt_1.5 (1553) minus
# colocalized_network_NPShr_gt_3 (184) = 1369, and dog_threshold_NPSd_gt_1.5 (1937) minus 184 =
# 1753, exactly. The panel itself doesn't display the hypergeometric p-value the caption
# mentions (that's `overlap_hypergeometric_p` in the same input file, computed but not part of
# this rendering) -- only the 3 counts.
#
# Previously undocumented as a real gap ("nothing renders the actual diagram from these
# numbers"), corrected once the real published panel turned out to be this simple -- no need for
# a 3-part diagram or matplotlib-venn after all.

suppressPackageStartupMessages({
  library(eulerr)
})

args <- commandArgs(trailingOnly = TRUE)
fig1e_counts_file <- args[1]  # fig1e_threshold_counts.tsv from CONSERVED_NETWORK
plot_out          <- args[2]

counts <- read.delim(fig1e_counts_file, header = FALSE, col.names = c("key", "value"))
get_val <- function(key) as.numeric(counts$value[counts$key == key])

n_human       <- get_val("human_threshold_NPSh_gt_1.5")
n_dog         <- get_val("dog_threshold_NPSd_gt_1.5")
n_colocalized <- get_val("colocalized_network_NPShr_gt_3")

n_human_only <- n_human - n_colocalized
n_dog_only   <- n_dog - n_colocalized

fit <- euler(c(
  "NPSh"     = n_human_only,
  "NPSd"     = n_dog_only,
  "NPSh&NPSd" = n_colocalized
))

p <- plot(fit,
  quantities = list(fontsize = 10),
  labels     = list(labels = c("NPSh > 1.5", "NPSd > 1.5"), fontsize = 9),
  fills      = list(fill = c("#3A89F9", "#D85A30", "#F26767"), alpha = 0.75),
  edges      = list(col = c("#3A89F9", "#D85A30"), lwd = 1.2)
)

pdf(plot_out, width = 5, height = 4)
print(p)
dev.off()
