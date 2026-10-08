#!/usr/bin/env Rscript
# plot_tissue_specificity.R -- formerly Fig. 4e (dropped from the figure 2026-10-06; the ORs are
# still reported in the Results text and Supplementary Table 13): odds ratio (Brain vs. other tissues) of dogCD and dogSize
# GWAS-region overlap with EpicDog cCREs.
#
# Source: testing_CRE_overlap_260617.R's "Tissue-specificity figure" section. That script builds a
# two-panel figure (raw overlap rates + this odds-ratio panel); the actual published Fig. 4e is
# just the odds-ratio panel alone (confirmed against the embedded figure in
# NATURE_Main_Figures.docx) -- the rates panel isn't part of the manuscript figure, so only the
# odds-ratio computation and plot are reproduced here.
#
# Uses epicdog_tissue_bp_overlaps.txt (real bedtools output, see check_epicdog_bp_overlap.sh),
# not the source script's original hardcoded literals -- those literals silently dropped
# mammary-gland tissue from the "other" pool (a mammary_gland/mammary filename mismatch in
# check_CRE_overlap.sh's original `cat` step; see data/Manuscript_draft/REVIEW_AND_SUGGESTIONS.md).
# Confirmed directly with KatarinaTe (2026-09-21): the manuscript's own published Fig. 4e panel
# still reflects the pre-fix numbers -- she has already produced a corrected version elsewhere and
# confirmed this script should use the fixed, live-recomputed values even though they don't match
# the currently-embedded manuscript image.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
epicdog_overlaps_file <- args[1]  # epicdog_tissue_bp_overlaps.txt: GWAS_set, Tissue_category, Overlap_bp, Total_bp
plot_out <- args[2]

epicdog_overlaps <- read.delim(epicdog_overlaps_file)

brain_cre_bp <- epicdog_overlaps$Total_bp[epicdog_overlaps$Tissue_category == "brain"][1]
other_cre_bp <- epicdog_overlaps$Total_bp[epicdog_overlaps$Tissue_category == "other"][1]

dogCD_O_brain <- epicdog_overlaps$Overlap_bp[epicdog_overlaps$GWAS_set == "dogCD" & epicdog_overlaps$Tissue_category == "brain"]
dogCD_O_other <- epicdog_overlaps$Overlap_bp[epicdog_overlaps$GWAS_set == "dogCD" & epicdog_overlaps$Tissue_category == "other"]
SIZE_O_brain  <- epicdog_overlaps$Overlap_bp[epicdog_overlaps$GWAS_set == "SIZE"  & epicdog_overlaps$Tissue_category == "brain"]
SIZE_O_other  <- epicdog_overlaps$Overlap_bp[epicdog_overlaps$GWAS_set == "SIZE"  & epicdog_overlaps$Tissue_category == "other"]

dogCD_tab <- matrix(c(dogCD_O_brain, brain_cre_bp - dogCD_O_brain,
                      dogCD_O_other, other_cre_bp - dogCD_O_other),
                    nrow = 2, byrow = TRUE,
                    dimnames = list(Tissue = c("Brain", "Other"), GWAS_Overlap = c("Yes", "No")))

SIZE_tab <- matrix(c(SIZE_O_brain, brain_cre_bp - SIZE_O_brain,
                     SIZE_O_other, other_cre_bp - SIZE_O_other),
                   nrow = 2, byrow = TRUE,
                   dimnames = list(Tissue = c("Brain", "Other"), GWAS_Overlap = c("Yes", "No")))

dogCD_fit <- fisher.test(dogCD_tab)  # two-sided, for a plottable finite 95% CI
SIZE_fit  <- fisher.test(SIZE_tab)

or_df <- data.frame(
  trait   = factor(c("dogCD", "dogSize"), levels = c("dogSize", "dogCD")),
  OR      = c(unname(dogCD_fit$estimate), unname(SIZE_fit$estimate)),
  CI_low  = c(dogCD_fit$conf.int[1], SIZE_fit$conf.int[1]),
  CI_high = c(dogCD_fit$conf.int[2], SIZE_fit$conf.int[2])
)

cat("dogCD OR/CI: ", dogCD_fit$estimate, dogCD_fit$conf.int, "\n")
cat("dogSize OR/CI: ", SIZE_fit$estimate, SIZE_fit$conf.int, "\n")

# Fixed x-axis range (0.8-1.5), shared with plot_region_forest.R (Fig. 4e), so the two
# odds-ratio forest plots sit on the same scale and can be shown/aligned together. KatarinaTe
# confirmed (2026-09-24) all real CI bounds fit inside this range. This replaces the per-plot
# dynamic range added 2026-09-23 (itself a fix for plot_region_forest.R silently clipping ACG's
# point off a stale hardcoded ceiling) -- that was a bug fix, not a rejection of fixed ranges in
# general. Warn (don't silently clip) if that stops being true, same failure mode as that bug.
x_limits <- c(0.8, 1.5)
if (any(or_df$CI_low < x_limits[1] | or_df$CI_high > x_limits[2])) {
  warning("plot_tissue_specificity.R: a trait's CI falls outside the fixed 0.8-1.5 x-axis range -- ",
          "it will be clipped in the plot. Check the real data before treating this figure as final.")
}

p <- ggplot(or_df, aes(y = trait)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey60", linewidth = 0.4) +
  geom_errorbar(aes(xmin = CI_low, xmax = CI_high, color = trait),
                orientation = "y", width = 0.25, linewidth = 1.0) +
  geom_point(aes(x = OR, color = trait), size = 2.5) +
  scale_color_manual(values = c("dogCD" = "#D55E00", "dogSize" = "grey40"), guide = "none") +
  scale_x_continuous(limits = x_limits) +
  labs(x = "Odds ratio (Brain vs. other tissues) Epic dog", y = NULL) +
  theme_classic(base_size = 12, base_family = "sans") +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())

ggsave(plot_out, p, width = 5.0, height = 1.8)
