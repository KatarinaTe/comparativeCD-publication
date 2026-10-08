#!/usr/bin/env Rscript
# ============================================================
# 02_plot_GO_clusters_all_cutoffs.R
#
# Plot the per-term tables written by 01_semantic_clustering.R (the
# top-community-filtered output), one figure per matching file.
#
# Input:
#   All files matching:
#   GO_with_without_dogCD.semantic_clustering.top_community.cut*.txt
#   (each has run1/run2 -log10 p, GO_semantic_cluster, cluster labels, and the
#    per-disease community columns produced by 01.)
#
# Output, for each input:
#   GO_run1_vs_run2_by_semantic_cluster.<suffix>.sigY.pdf   (suffix e.g. top_community.cut60)
#   GO_run1_vs_run2_by_semantic_cluster.<suffix>.sigY.png
# ============================================================

library(tidyverse)
library(cowplot)

#========================
# Settings
#========================

# -log10(p) threshold for calling a GO term significant in plots/counts
sig_cutoff <- 5

sig_suffix <- paste0("sig", gsub("\\.", "p", as.character(sig_cutoff)))

input_files <- list.files(
  path = ".",
  pattern = "^GO_with_without_dogCD\\.semantic_clustering\\.top_community\\..*\\.txt$",
  full.names = FALSE
)
if (length(input_files) == 0) {
  stop(
    "No input files found matching: ",
    "GO_with_without_dogCD.semantic_clustering.top_community.*.txt"
  )
}

message("Found ", length(input_files), " cutoff file(s):")
message(paste("  -", input_files, collapse = "\n"))

plot_one_cutoff_file <- function(input_file) {
  
  # capture everything between "semantic_clustering." and ".txt"
  # (e.g. "top_community.cut60")
  cluster_suffix <- stringr::str_match(
    basename(input_file),
    "semantic_clustering\\.(.*)\\.txt$"
  )[, 2]
  
  if (is.na(cluster_suffix) || cluster_suffix == "") {
    stop("Could not parse cutoff suffix from input file: ", input_file)
  }
  
  plot_outfile <- paste0(
    "GO_run1_vs_run2_by_semantic_cluster.",
    cluster_suffix,
    ".",
    sig_suffix,
    ".pdf"
  )
  
  message("\n========================================")
  message("Processing: ", input_file)
  message("Cluster suffix: ", cluster_suffix)
  message("Output PDF: ", plot_outfile)
  message("========================================")
  
  #========================
  # Read GO table
  #========================
  
  raw <- readr::read_tsv(input_file, show_col_types = FALSE)
  
  go_dat <- raw %>%
    mutate(
      run1 = as.numeric(run1),
      run2 = as.numeric(run2)
    ) %>% filter(run1 >= 0 | run2 >= 0)
  lost <- go_dat %>% ungroup() %>% filter(run1 >= sig_cutoff & run2 < sig_cutoff) %>% mutate(diff=run1-run2) %>% group_by(human_disease,GO_semantic_cluster,cluster_label) %>% summarize(nterms_lost=n(),max_diff=max(diff))
  hd <- go_dat %>% ungroup() %>% filter(run2 > sig_cutoff) %>% group_by(human_disease,GO_semantic_cluster,cluster_label) %>% summarize(nterms_hd=n())
  change <- hd %>% dplyr::full_join(lost, by = c("human_disease", "GO_semantic_cluster", "cluster_label")) %>% replace_na(list(nterms_lost=0,max_diff=0)) %>%
    mutate(frac=nterms_lost/(nterms_lost+nterms_hd)) %>% ungroup() %>% 
    dplyr::select(human=human_disease,nterms_hd,nterms_lost,frac,max_diff,cluster_label) %>% arrange(-frac,-nterms_lost)
  
  
  
  #========================
  # Select clusters to highlight
  #========================
  # Choose the 20 clusters represented by the highest run2 terms
  # after ranking terms within each disease.
  
  message(paste0("Clusters in ",input_file,": ",length(unique(go_dat$GO_semantic_cluster))))
  top_clusters <- go_dat %>% filter(run2 > sig_cutoff) %>% 
    group_by(human_disease,GO_semantic_cluster) %>% 
    summarize(run2=max(run2)) %>% ungroup() %>% 
    arrange(-run2) %>% group_by(human_disease) %>% mutate(rank=row_number()) %>% 
    ungroup() %>% mutate(rank=if_else(human_disease=="OCD",rank,rank+0.1)) %>% 
    group_by(GO_semantic_cluster) %>% summarize(rank=min(rank)) %>% 
    ungroup() %>% arrange(rank) %>%
    slice_head(n = 15) %>%
    mutate(
      cluster_num = as.numeric(str_extract(GO_semantic_cluster, "[0-9]+$"))
    ) %>% dplyr::select(-rank) %>% 
    arrange(cluster_num, GO_semantic_cluster)
  
  
  #========================
  # Build cluster-label table
  #========================
  
  all_clusters <- raw %>%
    dplyr::select(
      GO_semantic_cluster,
      cluster_label,
      cluster_label_wrapped,
      GO_id
    ) %>%
    distinct() %>%
    group_by(
      GO_semantic_cluster,
      cluster_label,
      cluster_label_wrapped
    ) %>%
    summarize(n_terms = n(), .groups = "drop")
  
  #========================
  # Define plotting colors and shapes
  #========================
  
  cluster_colors <- c(
    "#D55E00", "#0072B2", "#009E73", "#CC79A7", "#E69F00",
    "#56B4E9", "#8B4513", "#6A3D9A", "#B15928", "#1B9E77",
    "#E7298A", "#A6761D", "#66A61E", "#A50F15", "#2171B5",
    "#7570B3", "#D95F02", "#1F78B4", "#33A02C", "#B2182B",
    "#984EA3", "#4DAF4A", "#FF7F00", "#377EB8", "#E41A1C"
  )
  
  cluster_shapes <- c(
    17, 15, 18, 0, 1,
    2, 5, 6, 7, 8,
    9, 10, 12, 13, 14,
    17, 15, 18, 0, 1,
    2, 5, 6, 7, 8
  )
  
  if (nrow(top_clusters) > length(cluster_colors)) {
    stop("More highlighted clusters than available colors/shapes.")
  }
  
  # One key dataframe controls plot and legend aesthetics.
  cluster_key <- top_clusters %>%
    left_join(all_clusters, by = "GO_semantic_cluster") %>%
    mutate(
      cluster_label = if_else(
        is.na(cluster_label_wrapped) | cluster_label_wrapped == "",
        GO_semantic_cluster,
        cluster_label_wrapped
      ),
      legend_label = paste0(
        str_replace_all(cluster_label, "__", "\n")
      ),
      cluster_color = cluster_colors[seq_len(n())],
      cluster_shape = cluster_shapes[seq_len(n())]
    ) %>%
    dplyr::select(
      GO_semantic_cluster,
      cluster_num,
      cluster_label,
      legend_label,
      cluster_color,
      cluster_shape
    )
  
  legend_breaks <- cluster_key$GO_semantic_cluster
  legend_labels <- stats::setNames(cluster_key$legend_label, cluster_key$GO_semantic_cluster)
  cluster_color_map <- stats::setNames(cluster_key$cluster_color, cluster_key$GO_semantic_cluster)
  cluster_shape_map <- stats::setNames(cluster_key$cluster_shape, cluster_key$GO_semantic_cluster)
  
  stopifnot(identical(names(cluster_color_map), legend_breaks))
  stopifnot(identical(names(cluster_shape_map), legend_breaks))
  stopifnot(identical(names(legend_labels), legend_breaks))
  stopifnot(!anyNA(cluster_color_map))
  stopifnot(!anyNA(cluster_shape_map))
  stopifnot(!anyNA(legend_labels))
  
  n_legend <- length(legend_breaks)
  legend_ncol <- 5
  legend_height <- 0.4
  
  #========================
  # Add stable plotting ID
  #========================
  
  go_plot_dat <- go_dat %>%
    left_join(
      cluster_key %>%
        mutate(plot_cluster = GO_semantic_cluster) %>%
        dplyr::select(GO_semantic_cluster, plot_cluster),
      by = "GO_semantic_cluster"
    ) %>%
    mutate(
      plot_cluster = factor(plot_cluster, levels = legend_breaks)
    )
  
  #========================
  # Plot one disease
  #========================
  
  
  
  
  plot_one_disease <- function(in_disease) {
    
    disease_label <- c(
      OCD = "OCD",
      DEP = "Depression",
      SCH = "Schizophrenia"
    )
    count_lost <- raw %>% mutate(set=if_else(run1>5&run2>5,"both",if_else(run1>5 & run2 <=5,"lost",if_else(run1<=5 & run2 >5,"found","neither")))) %>% dplyr::select(human_disease,GO_id,GO_term,set) 
    
    
    counts <- raw  %>% filter(run1 >= sig_cutoff | run2 >= sig_cutoff) %>% filter(human_disease==in_disease) %>% 
      dplyr::select(GO_id,run1,run2) %>% distinct() %>% 
      pivot_longer(-GO_id) %>% mutate(sig=if_else(value >= sig_cutoff,TRUE,FALSE)) %>% 
      dplyr::select(-value) %>% pivot_wider(values_from=sig) %>% 
      group_by(run1,run2) %>% count()
    
    n_old <- sum((counts %>% filter(run1))$n) 
    n_new <- sum((counts %>% filter(!run1&run2))$n) 
    n_tot <- sum((counts %>% filter(run2))$n) 
    
    #filter(human_disease == in_disease&run1>5&run2 > sig_cutoff) %>% dplyr::select(GO_id) %>% distinct() %>% count()
    #n_tot <- go_plot_dat %>% filter(human_disease == in_disease&run2 > sig_cutoff) %>% dplyr::select(GO_id) %>% distinct() %>% count()
    
    n_new <- paste0("h+d: ",n_tot,"\n(new: ",n_new,")") 
    n_old <- paste0("h: ",n_old)
    
    pd_raw <- go_plot_dat %>%
      filter(human_disease == in_disease)

    maxX <- ceiling(max(pd_raw$run1, na.rm = TRUE) / 5) * 5
    maxY <- ceiling(max(pd_raw$run2, na.rm = TRUE) / 5) * 5
    jitter_x <- if (in_disease == "OCD") -maxX * 0.13 else -maxX * 0.08
    jitter_y <- -maxY * 0.05

    pd_all <- pd_raw %>%
      mutate(
        plot_run1 = if_else(run1 == 0, jitter_x, run1),
        plot_run2 = if_else(run2 == 0, jitter_y, run2)
      )

    pd_highlight <- pd_all %>%
      filter(!is.na(plot_cluster), run2 > sig_cutoff)

    pd_grey <- pd_all %>%
      filter(is.na(plot_cluster) | run2 <= sig_cutoff)

    # the data is filtered (in 01) to one top-level community per disease;
    # surface which community this panel represents.
    comm_lab <- pd_all %>%
      dplyr::filter(!is.na(community)) %>%
      dplyr::distinct(community) %>%
      dplyr::pull(community)
    comm_lab <- if (length(comm_lab)) paste(comm_lab, collapse = ", ") else NULL

    diag_max <- min(maxX, maxY)
    
    p <- ggplot(pd_all, aes(x = plot_run1, y = plot_run2))
    p <- p + geom_hline(yintercept = sig_cutoff, linetype = 2, linewidth = 0.25, color = "#8B2E2E")
    p <- p + geom_vline(xintercept = sig_cutoff, linetype = 2, linewidth = 0.25, color = "#8B2E2E")
    p <- p + geom_hline(yintercept = jitter_y, linetype = "dotted", linewidth = 0.25, color = "grey50")
    p <- p + geom_vline(xintercept = jitter_x, linetype = "dotted", linewidth = 0.25, color = "grey50")
    p <- p + annotate( "text", fontface="bold",x = maxX, y = sig_cutoff - (maxY / 100), label = n_old, hjust = 1, vjust = 1, size = 1.75, color = "grey20" ,lineheight=0.9 )
    p <- p + annotate( "text", fontface="bold",x = maxX, y = sig_cutoff + (maxY / 100), label = n_new, hjust = 1, vjust = 0, size = 1.75, color = "#8B2E2E" ,lineheight=0.9 )
    
    p <- p + annotate("segment", x = 0, y = 0, xend = diag_max, yend = diag_max, linewidth = 0.3, color = "grey50")
    p <- p + geom_point(data = pd_grey, color = "grey40", alpha = 0.25, shape = 16, size = 1)
    p <- p + geom_point(data = pd_highlight, aes(shape = plot_cluster, color = plot_cluster), alpha = 0.75, size = 2)
    
    p <- p + scale_color_manual(
      name = "Semantic clusters from GO terms (top 15)",
      values = cluster_color_map,
      breaks = legend_breaks,
      limits = legend_breaks,
      labels = legend_labels,
      drop = FALSE
    )
    
    p <- p + scale_shape_manual(
      name = "Semantic clusters from GO terms (top 15)",
      values = cluster_shape_map,
      breaks = legend_breaks,
      limits = legend_breaks,
      labels = legend_labels,
      drop = FALSE
    )
    
    p <- p + guides(
      color = guide_legend(
        ncol = legend_ncol,
        byrow = TRUE,
        title.position = "top",
        override.aes = list(
          shape = unname(cluster_shape_map[legend_breaks]),
          color = unname(cluster_color_map[legend_breaks]),
          size = 2,
          alpha = 0.75
        )
      ),
      shape = "none"
    )
    
    x_breaks <- pretty(c(0, maxX), n = 5)
    y_breaks <- pretty(c(0, maxY), n = 5)
    p <- p + scale_x_continuous(
      limits = c(jitter_x * 1.4, maxX),
      breaks = c(jitter_x, x_breaks),
      labels = c("NA", x_breaks)
    )
    p <- p + scale_y_continuous(
      limits = c(jitter_y * 1, maxY),
      breaks = c(jitter_y, y_breaks),
      labels = c("NA", y_breaks)
    )
    
    p <- p + labs(
      title = disease_label[[in_disease]],
     #subtitle = comm_lab,
      x = expression("human GWAS " * "(" * -log[10] * "p)"),
      y = expression("human + dogCD GWAS " * "(" * -log[10] * "p)")
    )
    
    p <- p + theme_cowplot(font_size = 7)
    p <- p + theme(
      plot.title = element_text(size=7,face = "bold"),
      plot.subtitle = element_text(size=6, color = "grey30"),
      legend.position = "none",
      axis.ticks = element_line(linewidth = 0.35),
      axis.ticks.length = unit(0.12, "cm"),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      axis.line = element_line(linewidth = 0.5)
    )
    
    p
  }
  
  #========================
  # Build complete legend from cluster_key
  #========================
  
  legend_df <- cluster_key %>%
    mutate(
      plot_cluster = factor(GO_semantic_cluster, levels = legend_breaks),
      x = seq_len(n()),
      y = 1
    )
  
  legend_plot <- ggplot(
    legend_df,
    aes(x = x, y = y, color = plot_cluster, shape = plot_cluster)
  )
  
  legend_plot <- legend_plot + geom_point(size = 2, alpha = 0.75)
  
  legend_plot <- legend_plot + scale_color_manual(
    name = "Semantic clusters from GO terms (top 15)",
    values = cluster_color_map,
    breaks = legend_breaks,
    limits = legend_breaks,
    labels = legend_labels,
    drop = FALSE
  )
  
  legend_plot <- legend_plot + scale_shape_manual(
    name = "Semantic clusters from GO terms (top 15)",
    values = cluster_shape_map,
    breaks = legend_breaks,
    limits = legend_breaks,
    labels = legend_labels,
    drop = FALSE
  )
  
  legend_plot <- legend_plot + guides(
    color = guide_legend(ncol = legend_ncol, byrow = TRUE, title.position = "top"),
    shape = guide_legend(ncol = legend_ncol, byrow = TRUE, title.position = "top")
  )
  
  legend_plot <- legend_plot + theme_void()
  legend_plot <- legend_plot + theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.box.just = "center",
    legend.title = element_text(size = 6.4, face = "bold"),
    legend.title.position = "left",
    legend.text = element_text(size = 5.9, lineheight = 0.9),
    legend.spacing.x = unit(0.07, "in"),
    legend.spacing.y = unit(0.03, "in"),
    legend.key.width = unit(0.163, "in"),
    legend.key.height = unit(0.112, "in"),
    legend.margin = margin(4, 4, 4, 4),
    legend.box.margin = margin(4, 4, 4, 4),
    legend.background = element_rect(color = "grey35", linewidth = 0.25, fill = "white")
  )
  
  legend <- cowplot::get_legend(legend_plot)
  
  #========================
  # Combine plots
  #========================
  
  diseases <- c("DEP", "SCH","OCD")
  
  plot_list <- purrr::map(diseases, plot_one_disease)
  
  for (i in 2:length(plot_list)) {
    plot_list[[i]] <- plot_list[[i]] + labs(y = NULL)
  }

  #========================
  # OCD NA-row detail strip
  #========================
  # Shows OCD points with run2 == 0 (no h+dogCD signal) spread along their
  # run1 x-axis value with y-jitter so overlapping points are distinguishable.

  # Points on the vertical NA line in OCD: run1 == 0 (no human-only signal).
  # Shown as a narrow strip to the LEFT of the OCD main panel, spread along
  # their run2 (y) values with x-jitter so overlapping points are visible.

  ocd_raw_strip      <- go_plot_dat %>% filter(human_disease == "OCD")
  ocd_maxY_strip     <- ceiling(max(ocd_raw_strip$run2, na.rm = TRUE) / 5) * 5
  ocd_jitter_y_strip <- -ocd_maxY_strip * 0.05   # matches jitter_y formula in plot_one_disease
  ocd_y_breaks_strip <- pretty(c(0, ocd_maxY_strip), n = 5)

  set.seed(42)
  strip_dat <- ocd_raw_strip %>%
    filter(run1 == 0) %>%
    mutate(
      plot_run2 = if_else(run2 == 0, ocd_jitter_y_strip, run2),
      x_jit     = runif(n(), -0.45, 0.45)
    )

  strip_grey      <- strip_dat %>% filter(is.na(plot_cluster))
  strip_highlight <- strip_dat %>% filter(!is.na(plot_cluster))

  ocd_strip <- ggplot(strip_dat, aes(x = x_jit, y = plot_run2)) +
    geom_hline(yintercept = sig_cutoff,         linetype = 2,        linewidth = 0.25, color = "#8B2E2E") +
    geom_hline(yintercept = ocd_jitter_y_strip, linetype = "dotted", linewidth = 0.25, color = "grey50") +
    geom_point(data = strip_grey,
               color = "grey40", alpha = 0.4, shape = 16, size = 1) +
    geom_point(data = strip_highlight,
               aes(color = plot_cluster, shape = plot_cluster),
               alpha = 0.75, size = 1.5) +
    scale_color_manual(values = cluster_color_map, breaks = legend_breaks,
                       limits = legend_breaks, drop = FALSE) +
    scale_shape_manual(values = cluster_shape_map, breaks = legend_breaks,
                       limits = legend_breaks, drop = FALSE) +
    scale_x_continuous(limits = c(-0.7, 0.7)) +
    scale_y_continuous(
      limits = c(ocd_jitter_y_strip * 1, ocd_maxY_strip),
      breaks = c(ocd_jitter_y_strip, ocd_y_breaks_strip),
      labels = c("NA", ocd_y_breaks_strip)
    ) +
    labs(x = "no h\nonly", y = NULL) +
    theme_cowplot(font_size = 7) +
    theme(
      legend.position  = "none",
      axis.text.y      = element_blank(),
      axis.ticks.y     = element_blank(),
      axis.line.y      = element_blank(),
      axis.text.x      = element_blank(),
      axis.ticks.x     = element_blank(),
      axis.title.x     = element_text(size = 5, color = "grey30", lineheight = 0.85),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      axis.line.x      = element_line(linewidth = 0.5)
    )

  maxwidth <- max(go_plot_dat$run1, na.rm = TRUE)

  relwidths <- go_plot_dat %>%
    group_by(human_disease) %>%
    summarize(width = max(run1, na.rm = TRUE), .groups = "drop") %>%
    mutate(
      width = width / maxwidth,
      width_inches = 5 * width,
      rel_width = width_inches / sum(width_inches),
      rel_width = if_else(human_disease == "OCD", rel_width + 0.2, rel_width)
    )

  ocd_rel_widths <- relwidths %>%
    filter(human_disease %in% diseases) %>%
    arrange(match(human_disease, diseases)) %>%
    pull(rel_width)

  # Nest the OCD strip (left) + OCD main (right) into one column,
  # aligned vertically so y-axes stay in register.
  ocd_col <- cowplot::plot_grid(
    ocd_strip,
    plot_list[[3]],
    nrow = 1,
    align = "h",
    axis = "tb",
    rel_widths = c(0.18, 0.82)
  )

  main_panel <- cowplot::plot_grid(
    plot_list[[1]],
    plot_list[[2]],
    ocd_col,
    nrow = 1,
    align = "h",
    axis = "tb",
    rel_widths = ocd_rel_widths
  )

  combined_p <- cowplot::plot_grid(
    main_panel,
    NULL,
    legend,
    ncol = 1,
    rel_heights = c(0.9, 0.01, legend_height)
  )
  # PDF
  ggsave(
    filename = plot_outfile,
    plot = combined_p,
    width = 6,
    height = 5.5,
    units = "in"
  )

  # PNG
  ggsave(
    filename = sub("\\.pdf$", ".png", plot_outfile),
    plot = combined_p,
    width = 6,
    height = 5.5,
    units = "in",
    dpi = 300
  )
  
  message("Wrote: ", plot_outfile)
  message("Wrote: ", sub("\\.pdf$", ".png", plot_outfile))
  
  invisible(combined_p)
}

plots <- purrr::map(input_files, plot_one_cutoff_file)

message("\nDone. Wrote plots for ", length(input_files), " cutoff file(s).")