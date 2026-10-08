#!/usr/bin/env Rscript
# ============================================================
# 03_OCD_NA_detail.R
#
# Supplementary plot: OCD GO terms that have h+dogCD signal
# but NO human-only signal (run1 = 0 / NA).
#
# Each row = one semantic cluster. Points = individual GO terms,
# jittered vertically within the row. X-axis = -log10p (h+dogCD).
# Colors and shapes match the main figure (02_plot_GO_clusters.v2.R).
#
# Input:
#   GO_with_without_dogCD.semantic_clustering.top_community.cut*.txt
# Output (one per input file):
#   OCD_NA_detail.<suffix>.pdf / .png
# ============================================================

library(tidyverse)
library(cowplot)

sig_cutoff <- 5

input_files <- list.files(
  path    = ".",
  pattern = "^GO_with_without_dogCD\\.semantic_clustering\\.top_community\\..*\\.txt$",
  full.names = FALSE
)
if (length(input_files) == 0) stop("No input files found.")
message("Found ", length(input_files), " file(s): ", paste(input_files, collapse = ", "))

# Identical color/shape vectors as the main figure
cluster_colors <- c(
  "#D55E00", "#0072B2", "#009E73", "#CC79A7", "#E69F00",
  "#56B4E9", "#8B4513", "#6A3D9A", "#B15928", "#1B9E77",
  "#E7298A", "#A6761D", "#66A61E", "#A50F15", "#2171B5",
  "#7570B3", "#D95F02", "#1F78B4", "#33A02C", "#B2182B",
  "#984EA3", "#4DAF4A", "#FF7F00", "#377EB8", "#E41A1C"
)
cluster_shapes <- c(
  17, 15, 18, 0, 1,
  2,  5,  6, 7, 8,
  9, 10, 12, 13, 14,
  17, 15, 18, 0, 1,
  2,  5,  6, 7, 8
)

plot_one_file <- function(input_file) {

  cluster_suffix <- stringr::str_match(
    basename(input_file), "semantic_clustering\\.(.*)\\.txt$"
  )[, 2]

  raw <- readr::read_tsv(input_file, show_col_types = FALSE)

  go_dat <- raw %>%
    mutate(run1 = as.numeric(run1), run2 = as.numeric(run2)) %>%
    filter(run1 >= 0 | run2 >= 0)

  # ---- Reproduce the same top-15 cluster selection as the main figure ----
  top_clusters <- go_dat %>%
    filter(run2 > sig_cutoff) %>%
    group_by(human_disease, GO_semantic_cluster) %>%
    summarize(run2 = max(run2), .groups = "drop") %>%
    arrange(-run2) %>%
    group_by(human_disease) %>%
    mutate(rank = row_number()) %>%
    ungroup() %>%
    mutate(rank = if_else(human_disease == "OCD", rank, rank + 0.1)) %>%
    group_by(GO_semantic_cluster) %>%
    summarize(rank = min(rank), .groups = "drop") %>%
    arrange(rank) %>%
    slice_head(n = 15) %>%
    mutate(cluster_num = as.numeric(str_extract(GO_semantic_cluster, "[0-9]+$"))) %>%
    dplyr::select(-rank) %>%
    arrange(cluster_num, GO_semantic_cluster)

  cluster_key <- top_clusters %>%
    left_join(
      raw %>% dplyr::select(GO_semantic_cluster, cluster_label) %>% distinct(),
      by = "GO_semantic_cluster"
    ) %>%
    mutate(
      display_label = if_else(is.na(cluster_label) | cluster_label == "",
                              GO_semantic_cluster, cluster_label),
      cluster_color = cluster_colors[seq_len(n())],
      cluster_shape = cluster_shapes[seq_len(n())]
    )

  legend_breaks     <- cluster_key$GO_semantic_cluster
  cluster_color_map <- setNames(cluster_key$cluster_color, cluster_key$GO_semantic_cluster)
  cluster_shape_map <- setNames(cluster_key$cluster_shape, cluster_key$GO_semantic_cluster)
  label_map         <- setNames(cluster_key$display_label,  cluster_key$GO_semantic_cluster)

  # ---- Filter to OCD, run1 == 0, run2 >= sig_cutoff ----
  ocd_na <- go_dat %>%
    filter(human_disease == "OCD", run1 == 0, run2 >= sig_cutoff) %>%
    left_join(
      cluster_key %>% dplyr::select(GO_semantic_cluster, display_label,
                                     cluster_color, cluster_shape),
      by = "GO_semantic_cluster"
    ) %>%
    mutate(
      in_top15     = GO_semantic_cluster %in% legend_breaks,
      row_label    = if_else(in_top15,
                             coalesce(display_label, cluster_label, GO_semantic_cluster),
                             "Other terms"),
      plot_cluster = if_else(in_top15, GO_semantic_cluster, NA_character_),
      plot_cluster = factor(plot_cluster, levels = legend_breaks)
    )

  if (nrow(ocd_na) == 0) {
    message("No OCD run1=0 rows with run2 >= ", sig_cutoff, " in ", input_file, ". Skipping.")
    return(invisible(NULL))
  }

  # Top-15 rows ordered by max run2 ascending (highest at top);
  # "Other terms" goes at the bottom.
  top15_order <- ocd_na %>%
    filter(in_top15) %>%
    group_by(row_label) %>%
    summarize(max_run2 = max(run2), .groups = "drop") %>%
    arrange(max_run2) %>%
    pull(row_label)

  ocd_na <- ocd_na %>%
    mutate(row_label = factor(row_label, levels = c("Other terms", top15_order)))

  n_rows      <- n_distinct(ocd_na$row_label)
  plot_height <- max(1.5, (n_rows * 0.28 + 1.2) / 2)

  p <- ggplot(ocd_na, aes(x = run2, y = row_label)) +
    geom_vline(xintercept = sig_cutoff, linetype = 2,
               linewidth = 0.3, color = "#8B2E2E") +
    # "Other terms" row: grey dots, no jitter
    geom_point(
      data  = ocd_na %>% filter(!in_top15),
      color = "grey50", alpha = 0.65, shape = 16, size = 1.8
    ) +
    # Top-15 clusters: colored + shaped, no jitter
    geom_point(
      data  = ocd_na %>% filter(in_top15),
      aes(color = plot_cluster, shape = plot_cluster),
      alpha = 0.85, size = 2
    ) +
    scale_color_manual(
      values = cluster_color_map, breaks = legend_breaks,
      limits = legend_breaks, drop = TRUE
    ) +
    scale_shape_manual(
      values = cluster_shape_map, breaks = legend_breaks,
      limits = legend_breaks, drop = TRUE
    ) +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.05))) +
    labs(
      x = expression("human + dogCD GWAS " * "(" * -log[10] * "p)"),
      y = NULL
    ) +
    theme_cowplot(font_size = 7) +
    theme(
      legend.position    = "none",
      axis.ticks         = element_line(linewidth = 0.35),
      panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.major.y = element_line(color = "grey92", linewidth = 0.2),
      panel.grid.minor   = element_blank(),
      axis.line          = element_line(linewidth = 0.4)
    )

  out_pdf <- paste0("OCD_NA_detail.", cluster_suffix, ".pdf")
  out_png <- sub("\\.pdf$", ".png", out_pdf)

  ggsave(out_pdf, plot = p, width = 6.5, height = plot_height, units = "in")
  ggsave(out_png, plot = p, width = 6.5, height = plot_height, units = "in", dpi = 300)
  message("Wrote: ", out_pdf, "  (", round(plot_height, 1), " in tall)")
}

purrr::walk(input_files, plot_one_file)
message("Done.")
