#!/usr/bin/env Rscript
# plot_item_factor_hist.R — factor-score histograms (row) + per-item response histograms (3 rows)
#
# Genericized via per-plot config tables (one row per plot) rather than 17 near-duplicated ggplot
# blocks — every functional/styling difference the original actually has (title, x-axis label,
# color, x-axis breaks/limits, title/axis-title margins, y-axis hjust, whether the y-axis is
# blanked) is preserved exactly per plot, including one genuine oddity reproduced as-is: item152
# has no `scale_x_continuous(breaks=1:6)` call unlike its 6 CCDF3 siblings.
#
# item147 uses `scale_x_continuous(limits = c(0.5, 6.5))` (a 2-element vector; the original's
# `c(0.5, 6, 5)` has 3 elements).
#
# N= annotations and the shared axis limits (x+y for the 3 factor-score plots, continuous mirt
# scores; y only for the 14 item plots, since their x-axis is the fixed 1-6 Likert instrument
# scale, kept literal) are computed from the real data here rather than hardcoded.

suppressMessages({
  library(ggplot2)
  library(ggtext)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
# args: f1_fam f2_fam f3_fam item7_fam item93_fam item95_fam item145_fam item146_fam item147_fam
#       item148_fam item149_fam item150_fam item151_fam item152_fam item153_fam item154_fam item155_fam
factor_fam_paths <- setNames(args[1:3], c("F1", "F2", "F3"))
item_ids <- c("7", "93", "95", "145", "146", "147", "148", "149", "150", "151", "152", "153", "154", "155")
item_fam_paths <- setNames(args[4:17], item_ids)

read_fam <- function(path, value_col) {
  d <- read.csv(path, sep = " ", header = FALSE, col.names = c("FID", "IID", "F", "M", "sex", value_col))
  d
}

factor_fams <- lapply(names(factor_fam_paths), function(fid) read_fam(factor_fam_paths[[fid]], fid))
names(factor_fams) <- names(factor_fam_paths)

item_fams <- lapply(item_ids, function(id) read_fam(item_fam_paths[[id]], paste0("item", id)))
names(item_fams) <- item_ids

# ── Factor-score histogram row (F1.hist + F2.hist + F3.hist) ───────────────────────────────────

factor_hist_configs <- list(
  list(id = "F1", color = "darkseagreen3", title = "CCD factor 1: Repetition severity",
       x_label = "", title_margin_left = 38, y_margin_left = 15, blank_y = FALSE),
  list(id = "F2", color = "coral1", title = "CCD factor 2: Compulsive staring/trancing/pacing",
       x_label = "", title_margin_left = 2, y_margin_left = 0, blank_y = TRUE),
  list(id = "F3", color = "dodgerblue2", title = "CCD factor 3: Repetition frequency",
       x_label = "Factor scores", title_margin_left = 38, y_margin_left = 15, blank_y = FALSE)
)

factor_binwidth <- 0.1
all_factor_scores  <- unlist(lapply(factor_fams, function(d) d[[6]]))
factor_range_pad   <- 0.1 * diff(range(all_factor_scores))
factor_shared_xlim <- range(all_factor_scores) + c(-factor_range_pad, factor_range_pad)
factor_breaks <- seq(factor_shared_xlim[1], factor_shared_xlim[2] + factor_binwidth, factor_binwidth)
factor_max_count <- max(sapply(factor_fams, function(d) {
  max(hist(d[[6]], breaks = factor_breaks, plot = FALSE)$counts)
}))
factor_shared_ylim <- c(0, ceiling(factor_max_count * 1.1))
factor_shared_breaks <- pretty(factor_shared_xlim)

make_factor_hist <- function(cfg) {
  d <- factor_fams[[cfg$id]]
  n_label <- paste0("N=", nrow(d))
  ggplot(d, aes(x = .data[[cfg$id]])) +
    geom_histogram(alpha = 2, binwidth = factor_binwidth, color = cfg$color, fill = cfg$color) +
    scale_fill_manual(values = c(cfg$color)) +
    scale_color_manual(values = c(cfg$color)) +
    ylim(factor_shared_ylim[1], factor_shared_ylim[2]) +
    scale_x_continuous(breaks = factor_shared_breaks, limits = factor_shared_xlim) +
    annotate("text", x = Inf, y = Inf, label = n_label, hjust = 1.3, vjust = 2.5, size = 4) +
    labs(title = cfg$title, x = cfg$x_label, y = "Number of dogs") +
    theme(plot.title.position = "plot",
          plot.title = element_textbox_simple(size = 10, padding = margin(5.5, 5.5, 5.5, 5.5),
                                               margin = margin(0, 0, 5.5, cfg$title_margin_left), fill = "cornsilk"),
          axis.title.x = element_textbox_simple(width = NULL, padding = margin(4, 4, 4, 4), margin = margin(4, 0, 0, 0),
                                                 fill = "white"),
          axis.title.y = element_textbox_simple(hjust = 0.65, orientation = "left-rotated", minwidth = unit(1, "in"), maxwidth = unit(2, "in"),
                                                 padding = margin(4, 4, 2, 4), margin = margin(0, 0, 2, cfg$y_margin_left), fill = "white")) +
    { if (cfg$blank_y) theme(axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank()) }
}

factor_plots <- lapply(factor_hist_configs, make_factor_hist)

# ── Per-item response histogram rows (CCDF1: 4 items, CCDF2: 3 items, CCDF3: 7 items) ──────────

item_hist_configs <- list(
  # CCDF1's items
  list(id = "7",   color = "darkseagreen3", title = "item7: Excitement can lead DOG to fixed repetitive behavior",
       x_label = "", x_breaks = NULL, x_limits = NULL, title_margin_left = 42.5, y_hjust = 1, y_margin_left = 20, blank_y = FALSE),
  list(id = "153", color = "darkslateblue",  title = "item153: It can be difficult to get DOG's attention when HE is engaged in a repetitive behavior",
       x_label = "", x_breaks = NULL, x_limits = NULL, title_margin_left = 0, y_hjust = 1, y_margin_left = 0, blank_y = TRUE),
  list(id = "154", color = "coral1",         title = "item154: It can be difficult to interrupt DOG when HE is engaged in a repetitive behavior",
       x_label = "Response", x_breaks = NULL, x_limits = NULL, title_margin_left = 42.5, y_hjust = 1, y_margin_left = 20, blank_y = FALSE),
  list(id = "155", color = "dodgerblue2",    title = "item155: One or more repetitive behaviors interfere with DOG's life",
       x_label = "Response", x_breaks = NULL, x_limits = NULL, title_margin_left = 0, y_hjust = 1, y_margin_left = 0, blank_y = TRUE),
  # CCDF2's items
  list(id = "93",  color = "turquoise3", title = "item93: DOG paces up and down, walks in circles and/or wanders with no direction or purpose",
       x_label = "Response", x_breaks = NULL, x_limits = NULL, title_margin_left = 42.5, y_hjust = 1.5, y_margin_left = 0, blank_y = FALSE),
  list(id = "95",  color = "thistle3",   title = "item95: DOG stares blankly at the walls or floor",
       x_label = "Response", x_breaks = NULL, x_limits = NULL, title_margin_left = 2, y_hjust = 1, y_margin_left = 0, blank_y = TRUE),
  list(id = "150", color = "dodgerblue4", title = "item150: How much time in a normal day does DOG spend showing trance-like behaviour?",
       x_label = "Response", x_breaks = 1:6, x_limits = NULL, title_margin_left = 2, y_hjust = 1, y_margin_left = 0, blank_y = TRUE),
  # CCDF3's items
  list(id = "145", color = "turquoise3", title = "item145: How much time in a normal day does DOG sped false digging (scratching at floor/carpet/etc)?",
       x_label = "", x_breaks = 1:6, x_limits = NULL, title_margin_left = 42.5, y_hjust = 3, y_margin_left = 5, blank_y = FALSE),
  list(id = "146", color = "thistle3",   title = "item146: How much time in a normal day does DOG spend chasing HIS tail, or circling?",
       x_label = "", x_breaks = 1:6, x_limits = NULL, title_margin_left = 2, y_hjust = 1, y_margin_left = 0, blank_y = TRUE),
  list(id = "147", color = "dodgerblue4", title = "item147: How much time in a normal day does DOG spend catching invisible flies?",
       x_label = "", x_breaks = 1:6, x_limits = c(0.5, 6.5), title_margin_left = 2, y_hjust = 1, y_margin_left = 0, blank_y = TRUE),
  list(id = "148", color = "darkseagreen3", title = "item148: How much time in a normal day does DOG spend licking, chewing, or sucking on themself?",
       x_label = "", x_breaks = 1:6, x_limits = NULL, title_margin_left = 42.5, y_hjust = 3, y_margin_left = 5, blank_y = FALSE),
  list(id = "149", color = "darkslateblue", title = "item149: How much time in a normal day does DOG spend licking, chewing, or sucking on toys, people, or other pets?",
       x_label = "Response", x_breaks = 1:6, x_limits = NULL, title_margin_left = 2, y_hjust = 1, y_margin_left = 0, blank_y = TRUE),
  list(id = "151", color = "coral1", title = "item151: How much time in a normal day does DOG spend licking floors?",
       x_label = "Response", x_breaks = 1:6, x_limits = NULL, title_margin_left = 2, y_hjust = 1, y_margin_left = 0, blank_y = TRUE),
  # item152 deliberately has no x_breaks — the original's own block omits scale_x_continuous(),
  # unlike its 6 CCDF3 siblings above (reproduced as a genuine source oddity, not "fixed").
  list(id = "152", color = "dodgerblue2", title = "item152: How much time in a normal day does DOG spend chasing lights or reflections?",
       x_label = "Response", x_breaks = NULL, x_limits = NULL, title_margin_left = 42.5, y_hjust = 3, y_margin_left = 5, blank_y = FALSE)
)

item_binwidth <- 1
item_max_count <- max(sapply(item_hist_configs, function(cfg) {
  d <- item_fams[[cfg$id]]
  col <- paste0("item", cfg$id)
  max(hist(d[[col]], breaks = seq(min(d[[col]], na.rm = TRUE) - 0.5, max(d[[col]], na.rm = TRUE) + 0.5, item_binwidth), plot = FALSE)$counts)
}))
item_shared_ylim <- c(0, ceiling(item_max_count * 1.1))

make_item_hist <- function(cfg) {
  d <- item_fams[[cfg$id]]
  col <- paste0("item", cfg$id)
  n_label <- paste0("N=", sum(!is.na(d[[col]])))

  p <- ggplot(d, aes(x = .data[[col]])) +
    geom_histogram(alpha = 1, binwidth = item_binwidth, color = cfg$color, fill = cfg$color) +
    scale_fill_manual(values = c(cfg$color)) +
    scale_color_manual(values = c(cfg$color)) +
    ylim(item_shared_ylim[1], item_shared_ylim[2]) +
    annotate("text", x = Inf, y = Inf, label = n_label, hjust = 1.3, vjust = 2.5, size = 4) +
    labs(title = cfg$title, x = cfg$x_label, y = "Number of dogs") +
    theme(plot.title.position = "plot",
          plot.title = element_textbox_simple(size = 10, padding = margin(5.5, 5.5, 5.5, 5.5),
                                               margin = margin(0, 0, 5.5, cfg$title_margin_left), fill = "cornsilk"),
          axis.title.x = element_textbox_simple(width = NULL, padding = margin(4, 4, 4, 4), margin = margin(4, 0, 0, 0),
                                                 fill = "white"),
          axis.title.y = element_textbox_simple(hjust = cfg$y_hjust, orientation = "left-rotated", minwidth = unit(1, "in"), maxwidth = unit(2, "in"),
                                                 padding = margin(4, 4, 2, 4), margin = margin(0, 0, 2, cfg$y_margin_left), fill = "white")) +
    { if (cfg$blank_y) theme(axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank()) }

  if (!is.null(cfg$x_breaks)) {
    if (!is.null(cfg$x_limits)) {
      p <- p + scale_x_continuous(breaks = cfg$x_breaks, limits = cfg$x_limits)
    } else {
      p <- p + scale_x_continuous(breaks = cfg$x_breaks)
    }
  }
  p
}

item_plots <- setNames(lapply(item_hist_configs, make_item_hist), sapply(item_hist_configs, `[[`, "id"))

# ── Diagnostics ──────────────────────────────────────────────────────────────────────────────

sink("item_factor_hist_diagnostics.log", split = TRUE)
print(summary(factor_fams$F1$F1, exclude = NULL))
print(summary(factor_fams$F2$F2, exclude = NULL))
print(summary(factor_fams$F3$F3, exclude = NULL))
table(item_fams[["7"]]$item7)
table(item_fams[["147"]]$item147, exclude = NULL)
sink()

# ── Combined plots ───────────────────────────────────────────────────────────────────────────

png("factor_hist_row.png", width = 12, height = 5, units = "in", res = 150)
print(factor_plots[[1]] + factor_plots[[2]] + factor_plots[[3]])
dev.off()

png("item_hist_row_f1.png", width = 14, height = 5, units = "in", res = 150)
print(item_plots[["7"]] + item_plots[["153"]] + item_plots[["154"]] + item_plots[["155"]])
dev.off()

png("item_hist_row_f2.png", width = 11, height = 5, units = "in", res = 150)
print(item_plots[["93"]] + item_plots[["95"]] + item_plots[["150"]])
dev.off()

png("item_hist_row_f3.png", width = 20, height = 5, units = "in", res = 150)
print(item_plots[["145"]] + item_plots[["146"]] + item_plots[["147"]] + item_plots[["148"]] +
      item_plots[["149"]] + item_plots[["151"]] + item_plots[["152"]])
dev.off()
