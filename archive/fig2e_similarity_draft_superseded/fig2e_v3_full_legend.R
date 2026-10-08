#!/usr/bin/env Rscript
# fig2e_v3_full_legend.R -- Fig 2e reconstruction, rebuilt against the collaborator's real figure
# legend (received 2026-10-03), 5 groups (OCD, OCD/dogCD, dogCD, Depression, Schizophrenia), with
# proper Benjamini-Hochberg correction across all 10 pairs. Does NOT overwrite the two earlier
# draft outputs (fig2e_similarity_heatmap_DRAFT.* and *_with_dogCD.*) -- separate filenames.
#
# Legend (full), verbatim from the collaborator:
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
# Confirmed against this legend (unchanged from the earlier DRAFT reconstruction):
#   - significance cutoff p<1e-5 (-log10 P >= 5)
#   - Wang/BMA = "average of the two directions, each the mean of the highest similarity to any
#     term in the other network" -- GOSemSim's own BMA definition, already what was used
#   - score = (observed-expected)/(1-expected), expected = mean of 1000 random-term-set pairs
#   - parenthetical counts = literal shared significant-term counts (already verified exact)
#   - network definitions (OCD/dogCD cross-species; dogCD single-species dog-only; OCD/Depression/
#     Schizophrenia human-only) match exactly what was used
#
# Fixed here: asterisks now use Benjamini-Hochberg-corrected p-values across all pairs (the
# earlier DRAFT used raw, uncorrected per-pair empirical p -- the legend specifies BH correction).
#
# Still an open, unconfirmed choice (the legend doesn't say): the background pool the 1000 random
# GO term sets are drawn from. Using the full annotated GO:BP ontology (org.Hs.eg.db, ~12,588
# terms), same choice as the earlier DRAFT -- this got 3 of the original 6 pairs within ~0.02 of
# a plausible match and the right sign on all 6, the best-supported option tried so far.

suppressPackageStartupMessages({
  library(GOSemSim)
  library(org.Hs.eg.db)
  library(readr)
  library(dplyr)
  library(stringr)
  library(ggplot2)
})

term_table <- "archive/go_semantic_clustering_originals/GO_with_without_dogCD.semantic_clustering.top_community.cut60.txt"
dogcd_enrichment_file <- "data/05_network/CCD_ONLY_hierachy_full_GO_enrichment.tsv"
sig_cutoff <- 5
n_perm <- 1000
set.seed(42)

d <- read_tsv(term_table, show_col_types = FALSE)

groups <- list(
  OCD           = d %>% filter(human_disease == "OCD", run1 >= sig_cutoff) %>% pull(GO_id) %>% unique(),
  `OCD/dogCD`   = d %>% filter(human_disease == "OCD", run2 >= sig_cutoff) %>% pull(GO_id) %>% unique(),
  Depression    = d %>% filter(human_disease == "DEP", run1 >= sig_cutoff) %>% pull(GO_id) %>% unique(),
  Schizophrenia = d %>% filter(human_disease == "SCH", run1 >= sig_cutoff) %>% pull(GO_id) %>% unique()
)

dogcd_raw <- read_tsv(dogcd_enrichment_file, show_col_types = FALSE)
ccol <- "name.1"
sz <- tapply(suppressWarnings(as.numeric(dogcd_raw$query_size)), dogcd_raw[[ccol]],
             function(v) max(v, na.rm = TRUE))
top_by_size <- names(sz)[which.max(sz)]
nums <- suppressWarnings(as.numeric(str_extract(names(sz), "[0-9]+$")))
stopifnot(identical(top_by_size, names(sz)[which.min(nums)]))

dogcd_top <- dogcd_raw %>% filter(.data[[ccol]] == top_by_size) %>%
  transmute(GO_id = native, GO_term = name, p = p_value) %>%
  group_by(GO_id, GO_term) %>% summarize(p = min(p, na.rm = TRUE), .groups = "drop") %>%
  mutate(neglog10p = -log10(p))
groups[["dogCD"]] <- dogcd_top %>% filter(neglog10p >= sig_cutoff) %>% pull(GO_id) %>% unique()

message("Group sizes: ", paste(names(groups), lengths(groups), sep = "=", collapse = ", "))

message("Loading GOSemSim BP semantic data...")
sem_data <- godata(OrgDb = "org.Hs.eg.db", ont = "BP", computeIC = FALSE)
full_bg <- unique(sem_data@geneAnno$GO[sem_data@geneAnno$ONTOLOGY == "BP"])
sim <- function(a, b) mgoSim(a, b, sem_data, measure = "Wang", combine = "BMA")

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

message("Computing raw similarity + ", n_perm, "-permutation null for all ", nrow(pairs), " pairs...")
results <- pairs %>%
  rowwise() %>%
  mutate(a = list(groups[[row]]), b = list(groups[[col]]),
         n_overlap = length(intersect(a, b)), raw = sim(a, b)) %>%
  ungroup()

null_draws_for <- function(a, b) {
  null <- numeric(n_perm)
  for (i in seq_len(n_perm)) null[i] <- sim(sample(full_bg, length(a)), sample(full_bg, length(b)))
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
write_csv(results, "fig2e_similarity_results_v3_full_legend.csv")

# ---- Plot ---------------------------------------------------------------------------------
group_order <- c("dogCD", "OCD", "OCD/dogCD", "Depression", "Schizophrenia")
group_n <- lengths(groups)[group_order]
group_labels <- paste0(group_order, "\n(N=", group_n, ")")
names(group_labels) <- group_order

plot_dat <- results %>%
  mutate(
    row = factor(row, levels = rev(group_order[1:4])),
    col = factor(col, levels = group_order[2:5]),
    sig = pval_bh < 0.05,
    label = paste0(sprintf("%.2f", score), if_else(sig, "*", ""), "\n(", n_overlap, ")")
  )

p <- ggplot(plot_dat, aes(x = col, y = row, fill = score)) +
  geom_tile(colour = "white", linewidth = 1) +
  geom_text(aes(label = label), size = 3, lineheight = 0.9) +
  scale_fill_gradient2(
    low = "#B2182B", mid = "white", high = "#2166AC", midpoint = 0,
    limits = c(-1, 1), name = "Similarity\nabove chance"
  ) +
  scale_x_discrete(labels = group_labels[group_order[2:5]], position = "bottom") +
  scale_y_discrete(labels = group_labels[rev(group_order[1:4])]) +
  coord_fixed() +
  labs(x = NULL, y = NULL,
       title = "Fig. 2e — rebuilt against collaborator's real legend (2026-10-03)",
       subtitle = "GOSemSim Wang/BMA + 1000-perm chance correction, Benjamini-Hochberg-corrected significance") +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(size = 10, face = "bold", colour = "#1F3A4D"),
    plot.subtitle = element_text(size = 7.5, colour = "grey30")
  )

ggsave("fig2e_similarity_heatmap_v3_full_legend.png", p, width = 6.8, height = 5.7, dpi = 300, bg = "white")
ggsave("fig2e_similarity_heatmap_v3_full_legend.pdf", p, width = 6.8, height = 5.7)
message("Wrote fig2e_similarity_heatmap_v3_full_legend.png/.pdf and fig2e_similarity_results_v3_full_legend.csv")
