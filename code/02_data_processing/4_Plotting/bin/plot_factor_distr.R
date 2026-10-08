#!/usr/bin/env Rscript
# plot_factor_distr.R — F1/F2/F3 factor-score distribution histograms, stacked vertically
#
# Genericized across F1/F2/F3 via a small per-factor table — the three original blocks are
# identical except: the title/x-axis label text, the N= annotation, and whether axis.title.x is
# actively styled (only F3's block has it — the bottom plot in the vertical stack is the only one
# that needs its x-axis title styled). Every other margin/padding value is identical across all
# three and reproduced exactly.
#
# N (via nrow() — each QC6 fam already excludes missing-phenotype dogs via plink's --prune, so
# nrow() is exactly the plotted N) and the shared xlim/ylim (via the real combined range/max bin
# height across all three factors, with a small pad, since one shared range is applied across all
# three plots for visual comparability) are computed from the real data rather than hardcoded.

suppressMessages({
  library(ggplot2)
  library(ggtext)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
f1_fam_path <- args[1]
f2_fam_path <- args[2]
f3_fam_path <- args[3]

read_fam <- function(path) {
  read.csv(path, sep = " ", header = FALSE, col.names = c("FID", "IID", "F", "M", "sex", "pheno"))
}

fams <- list(F1 = read_fam(f1_fam_path), F2 = read_fam(f2_fam_path), F3 = read_fam(f3_fam_path))

factor_configs <- list(
  list(factor_id = "F1", x_label = "F1 score", show_x_axis_title = FALSE),
  list(factor_id = "F2", x_label = "F2 score", show_x_axis_title = FALSE),
  list(factor_id = "F3", x_label = "F3 score", show_x_axis_title = TRUE)
)

# Shared axis ranges across all three factors, computed from the real data (see header comment).
all_pheno   <- unlist(lapply(fams, function(d) d$pheno))
range_pad   <- 0.1 * diff(range(all_pheno))
shared_xlim <- range(all_pheno) + c(-range_pad, range_pad)

binwidth   <- 0.1
breaks     <- seq(shared_xlim[1], shared_xlim[2] + binwidth, binwidth)
max_count  <- max(sapply(fams, function(d) max(hist(d$pheno, breaks = breaks, plot = FALSE)$counts)))
shared_ylim <- c(0, ceiling(max_count * 1.1))

make_plot <- function(cfg) {
  pheno <- fams[[cfg$factor_id]]$pheno
  n_label     <- paste0("N=", length(pheno))
  sd2_val     <- sd(pheno) * 2 + summary(pheno)[4]
  q3_val      <- summary(pheno)[5]
  mean_val    <- summary(pheno)[4]
  median_val  <- summary(pheno)[3]

  base_theme <- theme(
    plot.title.position = "plot",
    plot.title = element_textbox_simple(size = 10, padding = margin(5.5, 5.5, 5.5, 5.5), margin = margin(0, 0, 5.5, 38),
                                         fill = "cornsilk"),
    axis.title.y = element_textbox_simple(hjust = 0.65, orientation = "left-rotated", minwidth = unit(1, "in"), maxwidth = unit(2, "in"),
                                           padding = margin(4, 4, 2, 4), margin = margin(0, 0, 2, 15), fill = "white")
  )

  p <- ggplot(data.frame(pheno = pheno), aes(x = pheno)) +
    geom_histogram(color = "black", fill = "green", alpha = 0.5, binwidth = binwidth) +
    annotate("text", x = Inf, y = Inf, label = n_label, hjust = 1.3, vjust = 2.5, size = 4) +
    labs(title = paste0("CCD", cfg$factor_id), x = cfg$x_label, y = "Number of dogs") +
    xlim(shared_xlim[1], shared_xlim[2]) +
    ylim(shared_ylim[1], shared_ylim[2]) +
    geom_vline(aes(xintercept = sd2_val), color = "red", linetype = "dashed", size = 0.5) +
    geom_vline(aes(xintercept = mean_val), color = "purple", linetype = "dashed", size = 0.5) +
    geom_vline(aes(xintercept = median_val), color = "orange", linetype = "dashed", size = 0.5) +
    geom_vline(aes(xintercept = q3_val), color = "blue", linetype = "dashed", size = 0.5) +
    base_theme

  if (cfg$show_x_axis_title) {
    p <- p + theme(axis.title.x = element_textbox_simple(width = NULL, padding = margin(4, 4, 4, 4), margin = margin(4, 0, 0, 0),
                                                           fill = "white"))
  }
  p
}

sink("factor_distr_diagnostics.log", split = TRUE)
for (cfg in factor_configs) {
  pheno <- fams[[cfg$factor_id]]$pheno
  print((summary(pheno)[6]) / (sd(pheno)))
  print(sd(pheno))
}
sink()

plots <- lapply(factor_configs, make_plot)

png("factor_distr.png", width = 7, height = 10, units = "in", res = 150)
print(plots[[1]] / plots[[2]] / plots[[3]])
dev.off()
