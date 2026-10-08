#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tidyverse)
  library(qqman)
})

psig <- 4e-7
psug <- 1e-6

# Shared Manhattan + QQ plotting, factored out of two near-identical mlma/polmm blocks
# (originally hardcoded to SIZE_QC6_LOCO.loco.mlma and item153's POLMM output only).
plot_gwas_results <- function(results, label, out_prefix, ymax = NULL) {
  if (is.null(ymax)) ymax <- max(-log10(results$p), na.rm = TRUE) + 1
  breaks_seq <- pretty(c(0, ymax))

  ggdata <- results %>%
    filter(Chr < 39) %>%
    group_by(Chr) %>%
    summarise(chr_len = max(bp)) %>%
    mutate(tot = cumsum(as.numeric(chr_len)) - as.numeric(chr_len)) %>%
    select(-chr_len) %>%
    left_join(results, ., by = "Chr") %>%
    arrange(Chr, bp) %>%
    mutate(BPcum = bp + tot)

  axisdf <- ggdata %>%
    filter(Chr < 39) %>%
    group_by(Chr) %>%
    summarize(center = (max(BPcum) + min(BPcum)) / 2)

  manhattan_plot <- ggplot(ggdata, aes(x = BPcum, y = -log10(p))) +
    geom_point(aes(color = as.factor(Chr)), alpha = 0.7, size = 0.8) +
    scale_color_manual(values = rep(c("grey2", "darkgrey"), 39)) +
    geom_point(data = subset(ggdata, p < psig), color = "firebrick4", size = 0.85, alpha = 0.75) +
    geom_point(data = subset(ggdata, p < psug & p > psig), color = "dodgerblue4", size = 0.85, alpha = 0.75) +
    scale_x_continuous(expand = c(0, 0), label = axisdf$Chr, breaks = axisdf$center) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, ymax), breaks = breaks_seq, label = breaks_seq) +
    xlab(" ") +
    ylab("-log10(p)") +
    ggtitle(label) +
    geom_hline(yintercept = -log10(psig), linetype = "dashed", linewidth = 0.5, alpha = 0.5, color = "firebrick4") +
    geom_hline(yintercept = -log10(psug), linetype = "dashed", linewidth = 0.5, alpha = 0.5, color = "dodgerblue4") +
    theme_classic() +
    theme(
      plot.title = element_text(size = 12),
      legend.position = "none",
      text = element_text(size = 12),
      axis.text.x = element_text(size = rel(0.75)),
      axis.text.y = element_text(size = rel(0.75)),
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )

  ggsave(manhattan_plot, file = paste0(out_prefix, "_manhattan.jpg"), dpi = 300, width = 35, height = 10, units = "cm")

  z <- qnorm(results$p / 2)
  lambda <- round(median(z^2, na.rm = TRUE) / qchisq(0.5, df = 1), 3)

  eo <- tibble(
    e = -log10(ppoints(length(results$p))),
    o = -log10(sort(results$p, decreasing = FALSE))
  )

  qq_plot <- ggplot(eo, aes(x = e, y = o)) +
    geom_point(alpha = 0.5, size = 0.2) +
    geom_abline(col = "firebrick4") +
    geom_hline(yintercept = -log10(psig), linetype = "dashed", linewidth = 0.5, alpha = 0.5, color = "firebrick4") +
    geom_hline(yintercept = -log10(psug), linetype = "dashed", linewidth = 0.5, alpha = 0.5, color = "dodgerblue4") +
    scale_x_continuous(expand = c(0, 0), limits = c(0, ymax)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, ymax)) +
    xlab(expression(Expected ~ ~-log[10](p))) +
    ylab(expression(Observed ~ ~-log[10](p))) +
    ggtitle(paste0("Lambda=", lambda)) +
    theme_classic()

  ggsave(qq_plot, file = paste0(out_prefix, "_qq.jpg"), dpi = 300, width = 10, height = 10, units = "cm")

  message(label, ": lambda=", lambda, ", n_sig(p<", psig, ")=", sum(results$p < psig, na.rm = TRUE))
}

# ── MLMA-LOCO (04_gwas/1_mlma-loco) — one Manhattan/QQ pair per CCD factor + SIZE ──
# File naming matches 1_mlma-loco/process_definitions.nf's RUN_MLMA_LOCO output:
# "${qc6_bed.baseName}_LOCO.loco.mlma", e.g. CCDF1_DA_MERGED_GENCOVE_AXIOM_QC6_LOCO.loco.mlma
mlma_phenotypes <- c("CCDF1", "CCDF2", "CCDF3", "SIZE")

for (pheno in mlma_phenotypes) {
  mlma_file <- Sys.glob(paste0(pheno, "_*_QC6_LOCO.loco.mlma"))
  if (length(mlma_file) == 0) {
    warning("No mlma file found for ", pheno, " (expected ", pheno, "_*_QC6_LOCO.loco.mlma); skipping.")
    next
  }
  results <- read.table(mlma_file[1], header = TRUE)
  plot_gwas_results(results, label = paste0(pheno, "_loco"), out_prefix = paste0(pheno, "_loco"))
}

# ── POLMM (04_gwas/2_polmm) — one Manhattan/QQ pair per survey item ──
# File naming matches 2_polmm/process_definitions.nf's output:
# "modi_simuMarkerOutput_POLMM_item${item_id}_FULLGRM_fullGeno.txt"
# The 14 dogCD survey items — matches 2_polmm.nf's own getItemIds()
# (originally hardcoded to item153 only).
polmm_item_ids <- c(7, 93, 95, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155)

for (item_id in polmm_item_ids) {
  polmm_file <- paste0("modi_simuMarkerOutput_POLMM_item", item_id, "_FULLGRM_fullGeno.txt")
  if (!file.exists(polmm_file)) {
    warning("No polmm file found for item ", item_id, " (expected ", polmm_file, "); skipping.")
    next
  }
  results <- read.table(polmm_file, header = TRUE)
  results$p <- results$Pvalue
  plot_gwas_results(results, label = paste0("item", item_id, "_polmm"),
                     out_prefix = paste0("CCDitem", item_id, "_polmm"), ymax = 10)
}
