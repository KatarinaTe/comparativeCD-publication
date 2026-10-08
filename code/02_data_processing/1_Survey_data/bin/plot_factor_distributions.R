#!/usr/bin/env Rscript
# plot_factor_distributions.R — final per-factor score distribution plots
# Source: 02_data_processing/1_Survey_data/Darwin_dog-CD_3-factorsolution.R, lines 579-615
#
# Needs all three factors' CCD3F files together — a barrier after the per-factor fan-out.

suppressMessages({
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
f1_path <- args[1]
f2_path <- args[2]
f3_path <- args[3]

F1_CCD3F <- read.table(f1_path, header = TRUE, sep = "\t")
F2_CCD3F <- read.table(f2_path, header = TRUE, sep = "\t")
F3_CCD3F <- read.table(f3_path, header = TRUE, sep = "\t")

sink("distribution_diagnostics.log", split = TRUE)

(summary(F1_CCD3F$F1)[6]) / (sd(F1_CCD3F$F1))

sink()

SD2VARDE <- sd(F1_CCD3F$F1) * 2 + summary(F1_CCD3F$F1)[4]
Q3VARDE  <- summary(F1_CCD3F$F1)[5]
MEAN     <- summary(F1_CCD3F$F1)[4]
MEDIAN   <- summary(F1_CCD3F$F1)[3]

png("factor_distribution_F1.png", width = 7, height = 7, units = "in", res = 150)
ggplot(F1_CCD3F, aes(x = F1)) + geom_histogram(color = "black", fill = "green", alpha = 0.5, binwidth = 0.1) +
  geom_vline(aes(xintercept = SD2VARDE), color = "red", linetype = "dashed", size = 0.5) +
  geom_vline(aes(xintercept = MEAN), color = "purple", linetype = "dashed", size = 0.5) +
  geom_vline(aes(xintercept = MEDIAN), color = "orange", linetype = "dashed", size = 0.5) +
  geom_vline(aes(xintercept = Q3VARDE), color = "blue", linetype = "dashed", size = 0.5) + theme_minimal()
dev.off()

sink("distribution_diagnostics.log", append = TRUE, split = TRUE)

(summary(F2_CCD3F$F2)[6]) / (sd(F2_CCD3F$F2))

sink()

SD2VARDE <- sd(F2_CCD3F$F2) * 2 + summary(F2_CCD3F$F2)[4]
Q3VARDE  <- summary(F2_CCD3F$F2)[5]
MEAN     <- summary(F2_CCD3F$F2)[4]
MEDIAN   <- summary(F2_CCD3F$F2)[3]

png("factor_distribution_F2.png", width = 7, height = 7, units = "in", res = 150)
ggplot(F2_CCD3F, aes(x = F2)) + geom_histogram(color = "black", fill = "green", alpha = 0.5, binwidth = 0.1) +
  geom_vline(aes(xintercept = SD2VARDE), color = "red", linetype = "dashed", size = 0.5) +
  geom_vline(aes(xintercept = MEAN), color = "purple", linetype = "dashed", size = 0.5) +
  geom_vline(aes(xintercept = MEDIAN), color = "orange", linetype = "dashed", size = 0.5) +
  geom_vline(aes(xintercept = Q3VARDE), color = "blue", linetype = "dashed", size = 0.5) + theme_minimal()
dev.off()

sink("distribution_diagnostics.log", append = TRUE, split = TRUE)

(summary(F3_CCD3F$F3)[6]) / (sd(F3_CCD3F$F3))

sink()

SD2VARDE <- sd(F3_CCD3F$F3) * 2 + summary(F3_CCD3F$F3)[4]
Q3VARDE  <- summary(F3_CCD3F$F3)[5]
MEAN     <- summary(F3_CCD3F$F3)[4]
MEDIAN   <- summary(F3_CCD3F$F3)[3]

png("factor_distribution_F3.png", width = 7, height = 7, units = "in", res = 150)
ggplot(F3_CCD3F, aes(x = F3)) + geom_histogram(color = "black", fill = "green", alpha = 0.5, binwidth = 0.1) +
  geom_vline(aes(xintercept = SD2VARDE), color = "red", linetype = "dashed", size = 0.5) +
  geom_vline(aes(xintercept = MEAN), color = "purple", linetype = "dashed", size = 0.5) +
  geom_vline(aes(xintercept = MEDIAN), color = "orange", linetype = "dashed", size = 0.5) +
  geom_vline(aes(xintercept = Q3VARDE), color = "blue", linetype = "dashed", size = 0.5) + theme_minimal()
dev.off()
