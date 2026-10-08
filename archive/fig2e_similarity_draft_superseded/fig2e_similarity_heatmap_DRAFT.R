#!/usr/bin/env Rscript
# fig2e_similarity_heatmap_DRAFT.R -- DRAFT reconstruction of Fig. 2e (semantic similarity
# across enriched GO terms), built 2026-10-01 from the manuscript Results text alone, since the
# original script/notebook hasn't been located yet. NOT verified against the real source --
# see the caveats below and in the chat writeup before treating this as final.
#
# Manuscript text reproduced here for reference:
#   "Measuring semantic similarity across all enriched terms confirms this (Fig. 2e): the
#   OCD/dogCD network is highly similar to both the depression and schizophrenia networks (0.79
#   and 0.67 of the maximum possible similarity above chance; P = 0.001), more similar than
#   depression and schizophrenia are to each other (0.49). The human-only OCD network, by
#   contrast, is only modestly similar to depression (0.34; P = 0.001) and less similar than
#   random term sets to schizophrenia (-0.26) or to the OCD/dogCD network (-0.23)."
#
# ---- What's solid (deterministic, fully confirmed against the manuscript) --------------------
# Term sets (significant GO terms per group, from the archived
# archive/go_semantic_clustering_originals/GO_with_without_dogCD.semantic_clustering.
# top_community.cut60.txt, the same -log10(p) table 02_plot_GO_clusters.R already uses for
# Fig. 2c/d):
#   OCD        = OCD, run1 (human-only) >= 5   -> N=8   (matches panel's "OCD (N=8)")
#   OCD/dogCD  = OCD, run2 (human+dogCD) >= 5   -> N=75  (matches "OCD/dogCD (N=75)")
#   Depression = DEP, run1 (human-only) >= 5    -> N=129 (matches "Depression (N=129)")
#   Schizophrenia = SCH, run1 (human-only) >= 5 -> N=116 (matches "Schizophrenia (N=116)")
# The "(N)" shown under each cell's score in the published panel is the literal GO-term overlap
# count between that pair's two term sets -- verified to match exactly for all 6 pairs (0, 7, 0,
# 66, 46, 50).
#
# ---- What's reconstructed (NOT verified) ------------------------------------------------------
# The similarity measure: GOSemSim Wang/BMA (best-match average) semantic similarity between the
# two term sets -- this correctly reproduces the RANKING of all 6 reported pairs.
#
# The chance-correction: score = (raw - chance_mean) / (1 - chance_mean), where chance_mean comes
# from 1000 permutations drawing random term sets of the same sizes from the full annotated
# GO:BP universe (12,588 terms, via org.Hs.eg.db) -- NOT from the paper's own already-enriched
# term pool, which was tried first and gave worse/wrong-signed results (see chat writeup). This
# choice get 5/6 pairs within ~0.02-0.13 of the reported value, and the 3 strongest pairs land on
# empirical P = 0.0010 (1/1000) -- matching the manuscript's "P = 0.001" exactly, which is
# reassuring but not proof this is the real method.
#
# Known discrepancy: OCD-Depression computes to ~0.00 here (not significant), vs. the reported
# 0.34 (P = 0.001). The other 5 pairs are all correctly signed and reasonably close in magnitude.
# OCD's own 8 significant terms are all broad transcription-regulation GO terms (not
# neuro/synaptic, unlike the other 3 groups) -- plausibly inflating OCD's chance baseline in a
# way this reconstruction doesn't fully correct for. Flagged for KatarinaTe to sanity-check this
# one cell specifically if/when the real source script turns up.

suppressPackageStartupMessages({
  library(GOSemSim)
  library(org.Hs.eg.db)
  library(readr)
  library(dplyr)
  library(ggplot2)
})

term_table <- "archive/go_semantic_clustering_originals/GO_with_without_dogCD.semantic_clustering.top_community.cut60.txt"
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
stopifnot(identical(lengths(groups), c(OCD = 8L, `OCD/dogCD` = 75L, Depression = 129L, Schizophrenia = 116L)))

message("Loading GOSemSim BP semantic data...")
sem_data <- godata(OrgDb = "org.Hs.eg.db", ont = "BP", computeIC = FALSE)
full_bg <- unique(sem_data@geneAnno$GO[sem_data@geneAnno$ONTOLOGY == "BP"])

sim <- function(a, b) mgoSim(a, b, sem_data, measure = "Wang", combine = "BMA")

pairs <- tribble(
  ~row,         ~col,
  "OCD",        "OCD/dogCD",
  "OCD",        "Depression",
  "OCD",        "Schizophrenia",
  "OCD/dogCD",  "Depression",
  "OCD/dogCD",  "Schizophrenia",
  "Depression", "Schizophrenia"
)

message("Computing raw similarity + ", n_perm, "-permutation chance correction for ", nrow(pairs), " pairs...")
results <- pairs %>%
  rowwise() %>%
  mutate(
    a = list(groups[[row]]),
    b = list(groups[[col]]),
    n_overlap = length(intersect(a, b)),
    raw = sim(a, b)
  ) %>%
  ungroup()

null_mean <- function(a, b) {
  null <- numeric(n_perm)
  for (i in seq_len(n_perm)) {
    null[i] <- sim(sample(full_bg, length(a)), sample(full_bg, length(b)))
  }
  null
}

results <- results %>%
  rowwise() %>%
  mutate(
    null_draws = list(null_mean(a, b)),
    chance = mean(null_draws),
    score = (raw - chance) / (1 - chance),
    pval = (sum(null_draws >= raw) + 1) / (n_perm + 1)
  ) %>%
  ungroup() %>%
  select(row, col, n_overlap, raw, chance, score, pval)

print(results)
write_csv(results, "fig2e_similarity_results_DRAFT.csv")

# ---- Plot, matching the published panel's layout ----------------------------------------------
group_order <- c("OCD", "OCD/dogCD", "Depression", "Schizophrenia")
group_n <- lengths(groups)[group_order]
group_labels <- paste0(group_order, "\n(N=", group_n, ")")
names(group_labels) <- group_order

plot_dat <- results %>%
  mutate(
    row = factor(row, levels = rev(group_order[1:3])),
    col = factor(col, levels = group_order[2:4]),
    sig = pval < 0.05,
    label = paste0(sprintf("%.2f", score), if_else(sig, "*", ""), "\n(", n_overlap, ")")
  )

p <- ggplot(plot_dat, aes(x = col, y = row, fill = score)) +
  geom_tile(colour = "white", linewidth = 1) +
  geom_text(aes(label = label), size = 3, lineheight = 0.9) +
  scale_fill_gradient2(
    low = "#B2182B", mid = "white", high = "#2166AC", midpoint = 0,
    limits = c(-1, 1), name = "Similarity\nabove chance"
  ) +
  scale_x_discrete(labels = group_labels[group_order[2:4]], position = "bottom") +
  scale_y_discrete(labels = group_labels[rev(group_order[1:3])]) +
  coord_fixed() +
  labs(x = NULL, y = NULL,
       title = "Fig. 2e DRAFT reconstruction -- not verified against source",
       subtitle = "GOSemSim Wang/BMA + 1000-perm chance correction (full GO:BP background)") +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(size = 10, face = "bold", colour = "#8B2E2E"),
    plot.subtitle = element_text(size = 8, colour = "grey30")
  )

ggsave("fig2e_similarity_heatmap_DRAFT.png", p, width = 5.5, height = 4.5, dpi = 300, bg = "white")
ggsave("fig2e_similarity_heatmap_DRAFT.pdf", p, width = 5.5, height = 4.5)
message("Wrote fig2e_similarity_heatmap_DRAFT.png/.pdf and fig2e_similarity_results_DRAFT.csv")
