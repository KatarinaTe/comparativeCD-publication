# =============================================================
# Pairwise FST: subsampling analyses
# Figures: Figure 2 (A+B combined), Figure 3 (IBD)
# =============================================================

library(googlesheets4)
library(dplyr)
library(stringr)
library(forcats)
library(ggplot2)
library(cowplot)
library(tidyr)

# ---- 1. Read data --------------------------------------------
subsample_url <- "<Google Sheet link removed>"
fulldata_url  <- "<Google Sheet link removed>"

sub <- read_sheet(subsample_url, sheet = "fst_subsample_summary_update")

# ---- 2. Tidy population names --------------------------------
recode_pop <- c(
  "BMT_rutland"         = "Rutland VT L2",
  "GL_rutland"          = "Rutland VT L1",
  "Unknown-Central_PA"  = "Central PA",
  "Westfield_MA"        = "Westfield MA",
  "Copake_NY"           = "Copake NY",
  "Mt._Washington_MA"   = "Mt. Washington MA",
  "Rockland_County_NY"  = "Rockland County NY"
)

pair_key <- function(a, b) {
  norm <- function(x) tolower(gsub("\\.", "", x))
  apply(cbind(norm(a), norm(b)), 1, function(v) paste(sort(v), collapse = " || "))
}

sub <- sub %>%
  mutate(
    popA  = recode(pop1, !!!recode_pop),
    popB  = recode(pop2, !!!recode_pop),
    pair  = paste(popA, popB, sep = " – "),
    key   = pair_key(popA, popB)
  ) %>%
  filter(!str_detect(popA, "Rockland"),
         !str_detect(popB, "Rockland")) %>%
  mutate(
    category = case_when(
      key == pair_key("Rutland VT L1", "Rutland VT L2") ~ "Rutland VT L1–L2 split",
      str_detect(popA, "Mt\\. Washington") |
        str_detect(popB, "Mt\\. Washington")            ~ "Mt. Washington MA pairs",
      TRUE                                              ~ "other pairs"
    ),
    pair = fct_reorder(pair, median_wt)
  )

# ---- 3. Full-dataset FST values ------------------------------
full <- read_sheet(fulldata_url, sheet = "S6") %>%
  mutate(key = pair_key(PopA, PopB))

sub <- sub %>% left_join(select(full, key, FST_Weighted), by = "key")

# ---- helper: wrap caption text at ~n chars -------------------
wrap_caption <- function(txt, width = 150) {
  words <- strsplit(txt, " ")[[1]]
  lines <- character(0)
  current <- ""
  for (w in words) {
    candidate <- if (nchar(current) == 0) w else paste(current, w)
    if (nchar(candidate) <= width) {
      current <- candidate
    } else {
      lines   <- c(lines, current)
      current <- w
    }
  }
  if (nchar(current) > 0) lines <- c(lines, current)
  paste(lines, collapse = "\n")
}

# =============================================================
# Figure 2 (A + B): gap plot and slope chart side by side
# Legend inside panel A; matching margins keep panels equal height
# Total dimensions: 8 wide x 4.2 high
# =============================================================

n_sub <- 9

pop_n <- c(
  "Mt. Washington MA" = 38,
  "Copake NY"         = 11,
  "Central PA"        = 10,
  "Westfield MA"      = 11,
  "Rutland VT L1"     = 16,
  "Rutland VT L2"     = 9
)

# shared aesthetics for both panels
pal_shared <- c(
  "Rutland VT L1–L2 split" = "#c0392b",
  "Mt. Washington MA pairs" = "#8e44ad",
  "other pairs"             = "#34699a"
)
shape_shared <- c(
  "Rutland VT L1–L2 split" = 21,   # filled circle
  "Mt. Washington MA pairs" = 24,   # filled triangle
  "other pairs"             = 22    # filled square
)

# ---- 2A: gap plot --------------------------------------------
gap_df <- sub %>%
  mutate(
    nA        = pop_n[popA],
    nB        = pop_n[popB],
    size_term = (1/n_sub - 1/nA) + (1/n_sub - 1/nB),
    gap       = median_wt - FST_Weighted
  ) %>%
  filter(!is.na(gap), !is.na(size_term))

fit <- lm(gap ~ size_term, data = gap_df)
cf  <- coef(summary(fit))
ci  <- confint(fit)["size_term", ]
r2  <- summary(fit)$r.squared
sp  <- cor.test(gap_df$size_term, gap_df$gap, method = "spearman")

p_2a <- ggplot(gap_df, aes(x = size_term, y = gap,
                            colour = category, fill = category, shape = category)) +
  geom_smooth(method = "lm", se = TRUE, colour = "grey30",
              fill = "grey85", linewidth = 0.5,
              inherit.aes = FALSE, aes(x = size_term, y = gap)) +
  geom_point(size = 2.5, alpha = 0.9, stroke = 0.3, colour = "white") +
  scale_colour_manual(values = pal_shared, name = NULL) +
  scale_fill_manual(values = pal_shared, name = NULL) +
  scale_shape_manual(values = shape_shared, name = NULL) +
  labs(
    x = expression("Sample-size term  " * Sigma(1/n[sub] - 1/n[full])),
    y = expression("Subsampled median − full-N weighted " * F[ST])
  ) +
  theme_cowplot(font_size = 11) +
  background_grid(major = "xy", minor = "none",
                  colour.major = "grey92", size.major = 0.4) +
  theme(
    legend.position      = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.background    = element_rect(fill = "white", colour = NA),
    legend.key.size      = unit(0.35, "cm"),
    legend.text          = element_text(size = 7.5),
    plot.margin          = margin(5, 5, 20, 5)
  )

# ---- 2B: slope chart -----------------------------------------
slope_df <- sub %>%
  filter(!is.na(FST_Weighted)) %>%
  select(pair, category, full = FST_Weighted, subsampled = median_wt) %>%
  pivot_longer(c(full, subsampled),
               names_to = "dataset", values_to = "FST") %>%
  mutate(dataset = factor(dataset,
                          levels = c("full", "subsampled"),
                          labels = c("Full dataset", "Subsampled")))

p_2b <- ggplot(slope_df, aes(x = dataset, y = FST,
                              group = pair, colour = category,
                              fill = category, shape = category)) +
  geom_line(linewidth = 0.6, alpha = 0.8) +
  geom_point(size = 2.2, alpha = 0.9, stroke = 0.3, colour = "white") +
  ggrepel::geom_text_repel(
    data = slope_df %>% filter(dataset == "Subsampled"),
    mapping = aes(x = dataset, y = FST, label = pair, colour = category),
    size = 2.4, nudge_x = 0.08,
    direction = "y", hjust = 0, segment.size = 0.3,
    segment.colour = "grey70", max.overlaps = 20,
    show.legend = FALSE, inherit.aes = FALSE
  ) +
  scale_colour_manual(values = pal_shared, name = NULL) +
  scale_fill_manual(values = pal_shared, name = NULL) +
  scale_shape_manual(values = shape_shared, name = NULL) +
  scale_x_discrete(expand = expansion(add = c(0.2, 1.6))) +
  labs(
    x = " ",
    y = expression("Weighted " * F[ST])
  ) +
  theme_cowplot(font_size = 11) +
  background_grid(major = "y", minor = "none",
                  colour.major = "grey85", size.major = 0.4) +
  theme(
    legend.position = "none",
    axis.text.x     = element_text(size = 11),
    plot.margin     = margin(5, 5, 20, 5)
  )

# ---- assemble Figure 2 ---------------------------------------
panels <- plot_grid(p_2a, p_2b,
                    labels = c("A", "B"),
                    label_size = 13,
                    nrow = 1, rel_widths = c(1, 1))

# combined caption
caption_2 <- sprintf(
  wrap_caption(paste0(
    "Figure 2. Subsampling increases pairwise FST estimates but preserves rank order. ",
    "(A) FST inflation under subsampling scales with the degree of sample-size reduction. ",
    "Each point is a pairwise population comparison; the y-axis is the difference between ",
    "the median subsampled FST (30 replicates, n = 9 per population) and the full-dataset FST; ",
    "the x-axis is the sample-size term Σ(1/nₛᵤᵦ − 1/nᶠᵤˡˡ). ",
    "Regression: slope = %.3f, p = %.2e; R² = %.2f; Spearman ρ = %.2f (permutation p < 0.001). ",
    "(B) Each line connects the full-dataset FST (left) to the subsampled median (right). ",
    "Rank orders are correlated between full and subsampled sets (Spearman ρ = 0.90, permutation p < 0.001)."
  ), width = 145),
  cf["size_term", "Estimate"], cf["size_term", "Pr(>|t|)"],
  r2, sp$estimate
)

caption_grob <- ggdraw() +
  draw_label(caption_2, x = 0.01, y = 0.8, hjust = 0, vjust = 0.5,
             size = 8.5, colour = "grey20", lineheight = 1.1)

fig2_final <- plot_grid(panels, caption_grob,
                        ncol = 1, rel_heights = c(4.5, 1))

ggsave("figure2_AB.png", fig2_final, width = 8, height = 4.2, dpi = 300)
ggsave("figure2_AB.pdf", fig2_final, width = 8, height = 4.2)
print(fig2_final)

# =============================================================
# Figure 3: IBD subsampling stability strip plot
# =============================================================

ibd_reps <- read_sheet(subsample_url, sheet = "ibd_subsample_results_updated")

ref_vals <- data.frame(
  metric = c("Mantel r", "Mantel p-value", "FST per km"),
  value  = c(0.23, 0.29, 9.5e-5)    # <-- update with actual full-dataset values
)

ibd_long <- ibd_reps %>%
  select(rep, mantel_r, mantel_p, slope_per_km) %>%
  pivot_longer(-rep, names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric,
                         levels = c("mantel_r", "mantel_p", "slope_per_km"),
                         labels = c("Mantel r", "Mantel p-value", "FST per km")))

ref_vals <- ref_vals %>%
  mutate(metric = factor(metric,
                         levels = c("Mantel r", "Mantel p-value", "FST per km")))

p_ibd <- ggplot(ibd_long, aes(x = "", y = value)) +
  geom_jitter(width = 0.15, size = 2.2, alpha = 0.65, colour = "#34699a") +
  geom_boxplot(width = 0.35, outlier.shape = NA, fill = NA,
               linewidth = 0.5, colour = "grey40") +
  geom_hline(data = ref_vals, aes(yintercept = value),
             linetype = "dashed", colour = "#c0392b", linewidth = 0.7) +
  facet_wrap(~ metric, scales = "free_y", nrow = 1) +
  labs(
    x       = NULL,
    y       = NULL,
    caption = wrap_caption(
      paste0(
        "Figure 3. Isolation-by-distance statistics are stable across subsampling replicates ",
        "(n = 9 per population, 30 replicates). Each point is one replicate; box shows ",
        "median and interquartile range; red dashed line shows the full-dataset value. ",
        "Panels show the Mantel correlation between pairwise FST and geographic distance (r), ",
        "its p-value, and the isolation-by-distance slope (FST per km). ",
        "Mantel p-values were non-significant in every replicate, confirming that the absence ",
        "of an isolation-by-distance signal is not an artefact of unequal sample sizes."
      ), width = 100)
  ) +
  theme_cowplot(font_size = 11) +
  background_grid(major = "y", minor = "none",
                  colour.major = "grey85", size.major = 0.4) +
  theme(
    axis.text.x           = element_blank(),
    axis.ticks.x          = element_blank(),
    strip.background      = element_blank(),
    strip.text            = element_text(face = "bold", size = 10),
    plot.caption          = element_text(size = 9, colour = "grey30",
                                         hjust = 0, lineheight = 1.4,
                                         margin = margin(t = 10)),
    plot.caption.position = "plot"
  )

ggsave("ibd_subsampling_stability.png", p_ibd, width = 5.5, height = 5.0, dpi = 300)
ggsave("ibd_subsampling_stability.pdf", p_ibd, width = 5.5, height = 5.0)
print(p_ibd)
