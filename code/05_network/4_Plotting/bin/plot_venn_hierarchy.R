#!/usr/bin/env Rscript
# plot_venn_hierarchy.R -- Fig. 2b: area-proportional 3-circle Venn diagrams of in_hierarchy
# genes, comparing each of dogCD / OCD / OCD-dogCD against the human-only Depression and
# Schizophrenia hierarchies. Replaces the previous Panel B (stacked bar + Euler diagram, see
# archive/fig2b_euler_diagram_superseded/) -- the published figure now uses these proportional
# circle Venns instead.
#
# Source: code/07_additional_plots_and_analyses/Figure2_new_panels/venn_ST8.R (KatarinaTe's
# drop-in, 2026-10-01). Geometry/plotting logic unmodified; only the input is changed from a live
# Google Sheet read (Supplementary Table 8) to the equivalent deposited file,
# data/05_network/gene_level_network_basis_table.csv -- same underlying gene-level membership
# table, checked into the repo 2026-09-28 (see data/05_network/README.md). Verified the deposited
# file reproduces the exact same 7 region counts per comparison as the original script's own
# header comment and the published Fig. 2b image:
#   dogCD      432 517 226 74 16  8  5
#   OCD        399 511 290 71 49 14  8
#   OCD/dogCD  420 511 106 62 28 14 17
#
# The published figure only uses the OCD and OCD/dogCD Venns; the dogCD Venn is generated for
# completeness (matching the source script) but isn't part of the manuscript-assembled panel --
# same "extra output, manual assembly picks a subset" pattern as Fig. 2a/2b being separate files.
#
# Usage: plot_venn_hierarchy.R <gene_level_network_basis_table.csv>
# Output: venn_hierarchy_<set>.pdf / .png in the working directory, 40 mm tall, 300 dpi.

library(tidyverse)
library(cowplot)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: plot_venn_hierarchy.R <gene_level_network_basis_table.csv>")
gene_table_file <- args[1]

out_dir   <- "."
height_mm <- 40
dpi       <- 300

# ---- Read data ---------------------------------------------------------------
st8 <- read_csv(gene_table_file, col_types = cols(.default = "c")) %>%
  filter(!is.na(gene)) %>%
  mutate(across(-gene, ~ coalesce(toupper(str_trim(.x)) == "TRUE", FALSE)))

gene_set <- function(col) st8$gene[st8[[col]]]

dep <- gene_set("in_hierarchy_Depression")
scz <- gene_set("in_hierarchy_Schizophrenia")

comparisons <- tribble(
  ~column,                  ~label,      ~colour,
  "in_hierarchy_dogCD",     "dogCD",     "#d54040",
  "in_hierarchy_OCD",       "OCD",       "#0074e6",
  "in_hierarchy_OCD/dogCD", "OCD/dogCD", "#d04900"
)

# ---- Style -------------------------------------------------------------------
# Greys for Depression/Schizophrenia; the comparison set uses its own colour,
# a light tint of it where it overlaps one grey circle, and a darker shade
# where all three overlap.
region_fills <- function(colour) {
  c(
    "100" = "#8f8f8f",                                       # Depression only
    "010" = "#b3b3b3",                                       # Schizophrenia only
    "110" = "#e4e4e4",                                       # Depression + Schizophrenia
    "001" = colour,                                          # comparison set only
    "101" = colorRampPalette(c("white", colour))(100)[42],   # Depression + set
    "011" = colorRampPalette(c("white", colour))(100)[42],   # Schizophrenia + set
    "111" = colorRampPalette(c(colour, "black"))(100)[8]     # all three
  )
}
name_pt  <- 7    # set-name font size (pt)
count_pt <- 6    # count font size (pt)
small_pt <- 4.5  # count font size in very small regions (pt)

# ---- Geometry helpers --------------------------------------------------------
# Area of the lens where two circles (radii r1, r2, centres d apart) overlap
lens_area <- function(r1, r2, d) {
  if (d >= r1 + r2) return(0)
  if (d <= abs(r1 - r2)) return(pi * min(r1, r2)^2)
  a <- r1^2 * acos((d^2 + r1^2 - r2^2) / (2 * d * r1))
  b <- r2^2 * acos((d^2 + r2^2 - r1^2) / (2 * d * r2))
  a + b - 0.5 * sqrt((-d + r1 + r2) * (d + r1 - r2) * (d - r1 + r2) * (d + r1 + r2))
}

# Centre distance giving an overlap area equal to `ov`
circle_dist <- function(r1, r2, ov) {
  if (ov <= 0) return(1.02 * (r1 + r2))
  if (ov >= pi * min(r1, r2)^2) return(abs(r1 - r2))
  uniroot(function(d) lens_area(r1, r2, d) - ov,
          c(abs(r1 - r2) + 1e-9, r1 + r2 - 1e-9), tol = 1e-10)$root
}

# Outline of the intersection of the circles in `idx`, as a polygon.
# An intersection of circles is convex, so its outline is the convex hull of the
# circle-edge points that lie inside all the other circles in `idx`.
circle_intersection <- function(cx, cy, r, idx, n_pts = 2000) {
  theta <- seq(0, 2 * pi, length.out = n_pts + 1)[-1]
  pts <- map_dfr(idx, ~ tibble(x = cx[.x] + r[.x] * cos(theta),
                               y = cy[.x] + r[.x] * sin(theta)))
  inside_all <- reduce(idx, function(keep, i) {
    keep & ((pts$x - cx[i])^2 + (pts$y - cy[i])^2 <= r[i]^2 * (1 + 1e-9))
  }, .init = rep(TRUE, nrow(pts)))
  pts <- pts[inside_all, ]
  if (nrow(pts) < 3) return(NULL)
  pts[chull(pts$x, pts$y), ]
}

region_counts <- function(D, S, A) {
  tibble(
    region = c("100", "010", "001", "110", "101", "011", "111"),
    n = c(
      length(setdiff(D, union(S, A))),
      length(setdiff(S, union(D, A))),
      length(setdiff(A, union(D, S))),
      length(setdiff(intersect(D, S), A)),
      length(setdiff(intersect(D, A), S)),
      length(setdiff(intersect(S, A), D)),
      length(intersect(intersect(D, S), A))
    )
  )
}

# ---- Plot --------------------------------------------------------------------
make_venn <- function(column, label, colour, grid_n = 600) {
  A <- gene_set(column)
  r <- sqrt(c(length(dep), length(scz), length(A)) / pi)  # area = number of genes

  # Depression at the origin, comparison set to the right, Schizophrenia below
  d_DA <- circle_dist(r[1], r[3], length(intersect(dep, A)))
  d_DS <- circle_dist(r[1], r[2], length(intersect(dep, scz)))
  d_SA <- circle_dist(r[2], r[3], length(intersect(scz, A)))
  d_SA <- min(max(d_SA, abs(d_DA - d_DS) + 1e-6), d_DA + d_DS - 1e-6)  # keep a valid triangle
  sx <- (d_DA^2 + d_DS^2 - d_SA^2) / (2 * d_DA)
  sy <- -sqrt(max(d_DS^2 - sx^2, 0))
  cx <- c(0, sx, d_DA)
  cy <- c(0, sy, 0)

  pad  <- 0.03 * max(r)
  xlim <- c(min(cx - r), max(cx + r)) + c(-pad, pad)
  ylim <- c(min(cy - r), max(cy + r)) + c(-pad, pad)
  step <- max(diff(xlim), diff(ylim)) / grid_n

  # Vector shapes, painted in layers: whole circles first, then the pairwise
  # overlaps on top, then the three-way overlap. Each region ends up showing
  # the colour of the top layer that covers it.
  layers <- list("100" = 1, "010" = 2, "001" = 3,
                 "110" = c(1, 2), "101" = c(1, 3), "011" = c(2, 3),
                 "111" = c(1, 2, 3))
  shapes <- imap_dfr(layers, ~ {
    poly <- circle_intersection(cx, cy, r, .x)
    if (is.null(poly)) NULL else mutate(poly, region = .y)
  }) %>%
    mutate(region = factor(region, levels = names(layers)))  # sets draw order

  # Point grid used only to find label positions (not drawn): which circles
  # each point is in, and its distance to the nearest circle edge
  grid <-expand_grid(x = seq(xlim[1], xlim[2], by = step),
                     y = seq(ylim[1], ylim[2], by = step)) %>%
    mutate(
      dD = sqrt((x - cx[1])^2 + (y - cy[1])^2),
      dS = sqrt((x - cx[2])^2 + (y - cy[2])^2),
      dA = sqrt((x - cx[3])^2 + (y - cy[3])^2),
      region = paste0(as.integer(dD <= r[1]), as.integer(dS <= r[2]), as.integer(dA <= r[3])),
      edge_dist = pmin(abs(r[1] - dD), abs(r[2] - dS), abs(r[3] - dA))
    ) %>%
    filter(region != "000")

  counts <- region_counts(dep, scz, A)
  stopifnot(sum(counts$n[counts$region %in% c("001", "101", "011", "111")]) == length(A))

  # Label each region at its most interior point
  extent <- max(diff(xlim), diff(ylim))
  labels <- grid %>%
    group_by(region) %>%
    slice_max(edge_dist, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    inner_join(counts, by = "region") %>%
    filter(n > 0) %>%
    mutate(
      radius = case_when(region == "100" ~ r[1], region == "010" ~ r[2],
                         region == "001" ~ r[3], TRUE ~ NA_real_),
      name   = case_when(region == "100" ~ "Depression", region == "010" ~ "Schizophrenia",
                         region == "001" ~ label, TRUE ~ NA_character_),
      is_set = !is.na(name),
      count_y = if_else(is_set, y - 0.12 * radius, y),
      count_size = if_else(edge_dist < 0.025 * extent, small_pt, count_pt),
      name_y = y + 0.15 * radius,
      name_size = case_when(!is_set ~ NA_real_,
                            radius > 0.7 * max(r) ~ name_pt,
                            nchar(name) < 7 ~ name_pt - 1,
                            TRUE ~ name_pt - 2),
      # geom_text sizes are in mm, so convert from pt
      count_size = count_size / .pt,
      name_size  = name_size / .pt
    )

  missing <- setdiff(counts$region[counts$n > 0], labels$region)
  if (length(missing) > 0)
    warning(label, ": region(s) ", paste(missing, collapse = ", "),
            " have genes but are too small to draw")

  p <- ggplot() +
    geom_polygon(data = shapes, aes(x, y, group = region, fill = region), colour = NA) +
    geom_text(data = labels, aes(x, count_y, label = n, size = count_size)) +
    geom_text(data = filter(labels, is_set), aes(x, name_y, label = name, size = name_size)) +
    scale_fill_manual(values = region_fills(colour), guide = "none") +
    scale_size_identity() +
    coord_equal(xlim = xlim, ylim = ylim, expand = FALSE) +
    theme_cowplot(font_size = name_pt) +
    theme_nothing() +  # cowplot: drop axes, ticks, legend
    theme(plot.margin = margin(0, 0, 0, 0),
          plot.background = element_rect(fill = "white", colour = NA))

  width_mm <- height_mm * diff(xlim) / diff(ylim)
  stub <- file.path(out_dir, paste0("venn_hierarchy_", str_replace_all(label, "/", "_")))
  mm <- 1 / 25.4  # save_plot takes inches
  save_plot(paste0(stub, ".pdf"), p, device = "pdf",  # vector: polygons + text
            base_width = width_mm * mm, base_height = height_mm * mm)
  save_plot(paste0(stub, ".png"), p, dpi = dpi, bg = "white",
            base_width = width_mm * mm, base_height = height_mm * mm)

  message(label, ": ", paste(counts$region, counts$n, sep = "=", collapse = "  "))
  invisible(p)
}

pwalk(comparisons, make_venn)
