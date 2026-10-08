#!/usr/bin/env Rscript
# plot_region_forest.R -- Fig. 4e (Fig. 4f before 2026-10-06): odds ratio (dogCD vs. dogSize) of GWAS-region overlap with UU
# brain cCREs, across 8 dog brain regions, ACG highlighted.
#
# Source: testing_CRE_overlap_260617.R's final "forest_plot2_dogCD" section (the third of three
# near-duplicate rewrites in that script -- a three-panel figure adding a cCRE-size-per-region
# panel). The actual published panel (confirmed against the embedded figure in
# NATURE_Main_Figures.docx) is just the single odds-ratio panel -- so only that panel is
# reproduced here, not the size/rate panels the source script also built.
#
# dogCDbp/Size_bp (total bp of the dogCD/SIZE GWAS region sets) are computed directly from the
# GWAS region BED files themselves, not hardcoded literals as in the source script.
#
# No pooled/meta-analytic OR across regions -- just the 8 per-region Fisher's-test ORs already
# in the plot, written out as a table + plain-text range sentence (region_OR_summary.txt)
# for the manuscript Results text. See the per-region loop below for why a pooled figure was
# dropped.
#
# Region order and OR values (2026-09-24, real run against 06b-clump-regions-100kb's live output
# from historical GWAS data -- see REPRODUCIBILITY_AUDIT.md): ACG highest at OR ~1.47, thalamus/
# hypothalamus/occipital/temporal/striatum/frontal/cerebellum descending to ~1.18. Supersedes this
# header's earlier "ACG ~1.81" claim from 2026-09-23, which was a transcription error made when
# documenting an unverified in-progress run, not a real recomputation -- same cached
# BUILD_UU_REGION_BP_OVERLAPS inputs both times, so ~1.47 is the correct, reproducible value.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
dogcd_overlaps_file <- args[1]  # brain_region_bp_overlaps.txt: Region, Overlap_bp, Total_bp
size_overlaps_file  <- args[2]  # SIZE_brain_region_bp_overlaps.txt: same shape
dogcd_gwas_bp       <- as.numeric(args[3])  # total bp of the dogCD GWAS region set
size_gwas_bp        <- as.numeric(args[4])  # total bp of the SIZE GWAS region set
plot_out            <- args[5]
summary_out         <- args[6]  # per-region OR/CI table + a plain-text range sentence

dogcd <- read.delim(dogcd_overlaps_file)
size  <- read.delim(size_overlaps_file)

regions <- dogcd$Region
results <- data.frame(region = regions, OR = NA, CI_low = NA, CI_high = NA, p = NA)

for (i in seq_along(regions)) {
  mat <- matrix(c(dogcd$Overlap_bp[i], dogcd_gwas_bp - dogcd$Overlap_bp[i],
                  size$Overlap_bp[i],  size_gwas_bp  - size$Overlap_bp[i]),
                nrow = 2, byrow = TRUE,
                dimnames = list(c("dogCD", "size"), c("overlap", "no_overlap")))
  ft <- fisher.test(mat)
  results$OR[i]      <- unname(ft$estimate)
  results$CI_low[i]  <- ft$conf.int[1]
  results$CI_high[i] <- ft$conf.int[2]
  results$p[i]        <- ft$p.value
}
print(results)

# Per-region OR/CI table + a plain-text summary sentence (range + CIs), written for the
# manuscript Results text -- no pooled/meta-analytic OR here. A pooled figure was considered
# (Mantel-Haenszel with a shared dogCD/SIZE-total margin across regions, matching the original
# testing_CRE_overlap_260617.R script) but dropped: that script's MH design was itself an earlier,
# unrelated Claude suggestion, not an established/validated method, and KatarinaTe (2026-09-24)
# opted to just report the 8 per-region ORs already shown in the plot instead.
min_row <- results[which.min(results$OR), ]
max_row <- results[which.max(results$OR), ]
summary_text <- sprintf(
  "Odds ratios for dogCD vs. dogSize GWAS-region overlap with brain-region cCREs ranged from OR=%.2f (95%% CI %.2f-%.2f) in %s to OR=%.2f (95%% CI %.2f-%.2f) in %s, the highest-enrichment region.",
  min_row$OR, min_row$CI_low, min_row$CI_high, min_row$region,
  max_row$OR, max_row$CI_low, max_row$CI_high, max_row$region
)
cat(summary_text, "\n")
writeLines(c(
  "Per-region odds ratios (dogCD vs. dogSize GWAS-region overlap with brain-region cCREs):",
  capture.output(print(results, row.names = FALSE)),
  "",
  summary_text
), con = summary_out)

results$region <- factor(results$region, levels = results$region[order(results$OR)])
results$highlight <- ifelse(as.character(results$region) == "ACG", "ACG", "Other regions")

# Fixed x-axis range (0.8-1.5), shared with plot_tissue_specificity.R (formerly Fig. 4e), so the two
# odds-ratio forest plots sit on the same scale and can be shown/aligned together. KatarinaTe
# confirmed (2026-09-24) all real CI bounds fit inside this range. This replaces the per-plot
# dynamic range added 2026-09-23 -- that fix was for a real bug (a stale hardcoded c(0.8, 1.4)
# ceiling silently clipped ACG's point off-screen), not a rejection of fixed ranges in general.
# Warn (don't silently clip) if that stops being true, same failure mode as that earlier bug.
x_limits <- c(0.8, 1.5)
if (any(results$CI_low < x_limits[1] | results$CI_high > x_limits[2])) {
  warning("plot_region_forest.R: a region's CI falls outside the fixed 0.8-1.5 x-axis range -- ",
          "it will be clipped in the plot. Check the real data before treating this figure as final.")
}

p <- ggplot(results, aes(y = region)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey60", linewidth = 0.4) +
  geom_errorbar(aes(xmin = CI_low, xmax = CI_high, color = highlight),
                orientation = "y", width = 0.35, linewidth = 0.9) +
  geom_point(aes(x = OR, color = highlight), size = 2.0) +
  scale_color_manual(values = c("ACG" = "#D55E00", "Other regions" = "#2C5F7C"), guide = "none") +
  scale_x_continuous(limits = x_limits) +
  labs(x = "Odds ratio (dogCD vs. dogSize) UU", y = NULL) +
  theme_classic(base_size = 12, base_family = "sans") +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())

ggsave(plot_out, p, width = 5.0, height = 2.6)
