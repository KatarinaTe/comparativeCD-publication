#!/usr/bin/env Rscript
# depth_histograms.R — depth distribution diagnostic plots (all merged dogs vs. the final ALLFAM
# subset)
#
# `DA_merged` (fam_qc3modi joined with the depth file) is recomputed fresh here from the same two
# frozen inputs already used in Stage C (BUILD_EXCLUDE_LIST_ROUND1) — a deterministic join, so
# this gives an identical object without adding a new cross-stage output for one diagnostic plot.
#
# `DA_all$IID <- DA_all$sampleID` here has no `as.character()` cast: `fam_2584dogs.txt` (which
# ALLFAM is restricted to) contains only clean, all-numeric dogIDs, so both sides of this join
# infer as the same type.

suppressMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(ggtext)
})

args <- commandArgs(trailingOnly = TRUE)
fam_qc3modi_path <- args[1]
depth_path       <- args[2]
allfam_fam_path  <- args[3]

# Recreate Stage C's DA_merged (fam_qc3modi joined with depth, by IID)
fam <- read.csv(fam_qc3modi_path, sep = "\t", header = TRUE)
DA_all_for_merge <- read.csv(depth_path, sep = "\t", header = FALSE, col.names = c("sampleID", "meanDepthALL"))
DA_all_for_merge$IID <- as.character(DA_all_for_merge$sampleID)
DA_merged <- left_join(fam, DA_all_for_merge, by = "IID")

DA_all <- read.csv(depth_path, sep = "\t", header = FALSE, col.names = c("sampleID", "meanDepthALL"))
head(DA_all)
DA_all$IID <- DA_all$sampleID

sink("depth_histograms.log", split = TRUE)

summary(DA_all$meanDepthALL, exclude = NULL)

sink()

ALLFAM <- read.csv(allfam_fam_path, sep = " ", header = FALSE, col.names = c("FID", "IID", "F", "M", "sex", "phe"))

sink("depth_histograms.log", append = TRUE, split = TRUE)
head(ALLFAM)
sink()

depth_allfam <- left_join(ALLFAM, DA_all, by = "IID")

sink("depth_histograms.log", append = TRUE, split = TRUE)
summary(depth_allfam$meanDepthALL, exclude = NULL)
sd(!is.na(depth_allfam$meanDepthALL))
nrow(depth_allfam)
sink()

hist_1 <- ggplot(DA_merged, aes(x = meanDepthALL)) +
  geom_histogram(alpha = 1, binwidth = 0.01, color = "coral3", fill = "coral3") +
  scale_fill_manual(values = c("coral3")) +
  scale_color_manual(values = c("coral3")) +
  xlim(-0.1, 5) +
  annotate("text", x = Inf, y = Inf, label = "N=3285", hjust = 1.3, vjust = 2.5, size = 4) +
  labs(title = "ALL DARWIN depth all sites", x = "Depth", y = "Number of dogs") +
  theme(plot.title.position = "plot",
        plot.title = element_textbox_simple(size = 10, padding = margin(5.5, 5.5, 5.5, 5.5), margin = margin(0, 0, 5.5, 38),
                                             fill = "cornsilk"),
        axis.title.x = element_textbox_simple(width = NULL, padding = margin(4, 4, 4, 4), margin = margin(4, 0, 0, 0),
                                               fill = "white"),
        axis.title.y = element_textbox_simple(hjust = 0.65, orientation = "left-rotated", minwidth = unit(1, "in"), maxwidth = unit(2, "in"),
                                               padding = margin(4, 4, 2, 4), margin = margin(0, 0, 2, 15), fill = "white"))

hist_4 <- ggplot(depth_allfam, aes(x = meanDepthALL)) +
  geom_histogram(alpha = 1, binwidth = 0.01, color = "cornflowerblue", fill = "cornflowerblue") +
  scale_fill_manual(values = c("cornflowerblue")) +
  scale_color_manual(values = c("cornflowerblue")) +
  xlim(-0.1, 5) +
  annotate("text", x = Inf, y = Inf, label = "N=2211", hjust = 1.3, vjust = 2.5, size = 4) +
  labs(title = "ALLFAM DARWIN depth all sites, mean=0.83, sd=0.35", x = "Depth", y = "Number of dogs") +
  theme(plot.title.position = "plot",
        plot.title = element_textbox_simple(size = 10, padding = margin(5.5, 5.5, 5.5, 5.5), margin = margin(0, 0, 5.5, 38),
                                             fill = "cornsilk"),
        axis.title.x = element_textbox_simple(width = NULL, padding = margin(4, 4, 4, 4), margin = margin(4, 0, 0, 0),
                                               fill = "white"),
        axis.title.y = element_textbox_simple(hjust = 0.65, orientation = "left-rotated", minwidth = unit(1, "in"), maxwidth = unit(2, "in"),
                                               padding = margin(4, 4, 2, 4), margin = margin(0, 0, 2, 15), fill = "white"))

png("depth_histograms.png", width = 6, height = 8, units = "in", res = 150)
print(hist_1 / hist_4)
dev.off()
