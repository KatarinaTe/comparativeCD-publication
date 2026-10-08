#!/usr/bin/env Rscript
# plot_semantic_similarity.R -- Fig. 2e: semantic similarity between networks' significantly
# enriched GO Biological Process terms, relative to chance.
#
# Method confirmed against the collaborator's real figure legend and her own intermediate
# computed values (2026-10-03) -- see this directory's README.md Notes for the full
# verification: all 15 pairwise raw-similarity values among DEP/OCD/SCH run1/run2 matched her
# own numbers to 3 decimals, and the chance-corrected score for the OCD-Depression cell matched
# within permutation noise.
#
# Legend (collaborator's wording):
#   "Semantic similarity between the GO Biological Process terms significantly enriched
#   (p < 1x10^-5) in each pair of networks, relative to chance. Similarity is the average of the
#   two directions, each the mean over one network's terms of the highest Wang semantic similarity
#   to any term in the other network, so that both networks count equally regardless of size.
#   Values give the fraction of the possible above-chance similarity achieved,
#   (observed - expected) / (1 - expected), where expected is the mean for 1000 pairs of random GO
#   term sets of the same sizes (0, as similar as random terms; 1, identical sets; negative, less
#   similar than random). Axis labels give each network's total number of significant GO terms;
#   numbers in parentheses within cells give the number of significant GO terms shared by the two
#   networks. *, more similar than random (Benjamini-Hochberg-corrected empirical p < 0.05).
#   OCD/dogCD, colocalized cross-species network for OCD and dogCD; dogCD, single-species, dog-only
#   network; OCD, depression and schizophrenia, human-only networks."
#
# A note on why this ISN'T GOSemSim's own `combine="BMA"`: that option pools every term from both
# sets together before averaging, so it's size-weighted -- dominated by whichever network has more
# terms. The real method averages the two *directional* means with equal weight regardless of
# size (GOSemSim's `combine=NULL` full matrix, row-maxed and column-maxed separately, then
# averaged) -- only visibly different for badly size-imbalanced pairs like OCD (8 terms) vs.
# Depression (129 terms), which is exactly the pair every earlier draft reconstruction got wrong.
# The permutation background (what pool the random term sets are drawn from) is NOT specified by
# the legend; the full genome-wide GO:BP annotation (org.Hs.eg.db) is used here, confirmed to
# reproduce the collaborator's own chance/expected values within normal permutation noise.
#
# Usage: plot_semantic_similarity.R <per_term_table.tsv> <out_prefix>
#   per_term_table.tsv: SEMANTIC_CLUSTERING's GO_with_without_dogCD.semantic_clustering.
#                        top_community.cut60.txt (human_disease, run1, run2, GO_id columns)

suppressPackageStartupMessages({
  library(GOSemSim)
  library(org.Hs.eg.db)
  library(readr)
  library(dplyr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("Usage: plot_semantic_similarity.R <per_term_table.tsv> <out_prefix>")
per_term_file <- args[1]
out_prefix    <- args[2]

sig_cutoff <- 5       # -log10(P) >= 5, i.e. P < 1e-5
n_perm     <- 1000
set.seed(42)

d <- read_tsv(per_term_file, show_col_types = FALSE)

groups <- list(
  OCD           = d %>% filter(human_disease == "OCD", run1 >= sig_cutoff) %>% pull(GO_id) %>% unique(),
  `OCD/dogCD`   = d %>% filter(human_disease == "OCD", run2 >= sig_cutoff) %>% pull(GO_id) %>% unique(),
  dogCD         = d %>% filter(human_disease == "dogCD", run1 >= sig_cutoff) %>% pull(GO_id) %>% unique(),
  Depression    = d %>% filter(human_disease == "DEP", run1 >= sig_cutoff) %>% pull(GO_id) %>% unique(),
  Schizophrenia = d %>% filter(human_disease == "SCH", run1 >= sig_cutoff) %>% pull(GO_id) %>% unique()
)
message("Group sizes: ", paste(names(groups), lengths(groups), sep = "=", collapse = ", "))

message("Loading GOSemSim BP semantic data...")
sem_data <- godata(OrgDb = "org.Hs.eg.db", ont = "BP", computeIC = FALSE)
full_bg  <- unique(sem_data@geneAnno$GO[sem_data@geneAnno$ONTOLOGY == "BP"])

# Equal-weighted average of the two directional best-match means -- see header note on why this
# is NOT GOSemSim's own combine="BMA".
equal_weighted_similarity <- function(a, b) {
  m <- mgoSim(a, b, sem_data, measure = "Wang", combine = NULL)
  sim_a_to_b <- mean(apply(m, 1, max, na.rm = TRUE))
  sim_b_to_a <- mean(apply(m, 2, max, na.rm = TRUE))
  (sim_a_to_b + sim_b_to_a) / 2
}

pairs <- tribble(
  ~row,         ~col,
  "dogCD",      "OCD",
  "dogCD",      "OCD/dogCD",
  "dogCD",      "Depression",
  "dogCD",      "Schizophrenia",
  "OCD",        "OCD/dogCD",
  "OCD",        "Depression",
  "OCD",        "Schizophrenia",
  "OCD/dogCD",  "Depression",
  "OCD/dogCD",  "Schizophrenia",
  "Depression", "Schizophrenia"
)

message("Computing similarity + ", n_perm, "-permutation null for ", nrow(pairs), " pairs...")
results <- pairs %>%
  rowwise() %>%
  mutate(a = list(groups[[row]]), b = list(groups[[col]]),
         n_overlap = length(intersect(a, b)),
         raw = equal_weighted_similarity(a, b)) %>%
  ungroup()

null_draws_for <- function(a, b) {
  null <- numeric(n_perm)
  for (i in seq_len(n_perm)) {
    null[i] <- equal_weighted_similarity(sample(full_bg, length(a)), sample(full_bg, length(b)))
  }
  null
}

results <- results %>%
  rowwise() %>%
  mutate(
    null_draws = list(null_draws_for(a, b)),
    chance = mean(null_draws),
    score = (raw - chance) / (1 - chance),
    pval_raw = (sum(null_draws >= raw) + 1) / (n_perm + 1)
  ) %>%
  ungroup() %>%
  mutate(pval_bh = p.adjust(pval_raw, method = "BH")) %>%
  select(row, col, n_overlap, raw, chance, score, pval_raw, pval_bh)

print(results)
write_csv(results, paste0(out_prefix, "_results.csv"))

# ---- Plot -----------------------------------------------------------------------------------
# Layout and styling match panel e of the assembled Figure 2 (Figure2_261003.pdf, measured from its
# PDF 2026-10-06): rows dogCD / OCD / Depression / Schizophrenia, columns OCD / Depression /
# Schizophrenia / OCD/dogCD (upper triangle); 5pt Helvetica cell and axis text; 4pt legend; tiles
# 18.2 x 16.7pt; and a muted diverging palette whose end colours were fitted (in Lab, as ggplot
# interpolates) to that panel's 10 tile colours -- the red end rests on only two light cells, so
# its far end is the least certain. Only placement changes here: `pairs` above (and so the
# permutation draws) is untouched, each pair is just put in the upper triangle of this layout.
row_order <- c("dogCD", "OCD", "Depression", "Schizophrenia")
col_order <- c("OCD", "Depression", "Schizophrenia", "OCD/dogCD")
layout_order <- c("dogCD", "OCD", "Depression", "Schizophrenia", "OCD/dogCD")
group_n <- lengths(groups)
row_labels <- setNames(paste0(row_order, "\n(", group_n[row_order], ")"), row_order)
col_wrap   <- c(OCD = "OCD", Depression = "Depres-\nsion", Schizophrenia = "Schizo-\nphrenia",
                `OCD/dogCD` = "OCD/\ndogCD")
col_labels <- setNames(paste0(col_wrap[col_order], "\n(", group_n[col_order], ")"), col_order)

plot_dat <- results %>%
  mutate(
    first  = match(row, layout_order) < match(col, layout_order),
    r      = if_else(first, row, col),
    c      = if_else(first, col, row),
    r      = factor(r, levels = rev(row_order)),
    c      = factor(c, levels = col_order),
    sig    = pval_bh < 0.05,
    # a plain "-": the base pdf device's Helvetica encoding draws it as a true minus sign
    label  = paste0(sprintf("%.2f", score), if_else(sig, "*", ""), "\n(", n_overlap, ")"),
    text_colour = if_else(score > 0.6, "white", "grey15")
  )

p <- ggplot(plot_dat, aes(x = c, y = r, fill = score)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = label, colour = text_colour), size = 5 / .pt, lineheight = 0.85,
            family = "Helvetica") +
  scale_colour_identity() +
  scale_fill_gradient2(
    low = "#C53F3D", mid = "#F0EFEB", high = "#0B356B", midpoint = 0,
    limits = c(-1, 1), name = "Similarity\nabove chance",
    breaks = c(-1, 0, 1),
    guide = guide_colorbar(
      direction = "horizontal", title.position = "left",
      barwidth = unit(0.5, "in"), barheight = unit(0.06, "in"),
      title.hjust = 0, label.position = "bottom", order = 1
    )
  ) +
  scale_x_discrete(labels = col_labels, position = "bottom") +
  scale_y_discrete(labels = row_labels) +
  coord_fixed(ratio = 16.73 / 18.19, clip = "off") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 5, base_family = "Helvetica") +
  theme(
    panel.grid = element_blank(),
    axis.ticks.length.y = unit(0, "pt"),
    axis.text.x = element_text(size = 5, colour = "black", lineheight = 0.85, margin = margin(t = 1)),
    axis.text.y = element_text(size = 5, colour = "black", lineheight = 0.85, hjust = 1,
                               margin = margin(r = 1)),
    legend.title = element_text(size = 4, lineheight = 0.85),
    legend.text = element_text(size = 4),
    legend.position = "bottom",
    legend.box.margin = margin(t = -4),
    legend.margin = margin(0, 0, 0, 0),
    plot.margin = margin(1, 1, 1, 1)
  )

# Canvas sized to panel e's slot in Figure 2, so embedding needs ~1:1 scaling and these absolute
# pt font sizes survive unchanged.
ggsave(paste0(out_prefix, ".pdf"), p, width = 112 / 72.27, height = 104 / 72.27, device = grDevices::pdf,
       useDingbats = FALSE)  # base pdf device: real Helvetica metrics (not in this container as a system font)
ggsave(paste0(out_prefix, ".png"), p, width = 112 / 72.27, height = 104 / 72.27, dpi = 600, bg = "white")
message("Wrote ", out_prefix, ".{pdf,png} and ", out_prefix, "_results.csv")
