#!/usr/bin/env Rscript
# fig2e_add_dogcd.R -- extends the Fig 2e DRAFT reconstruction with a 5th group, dogCD
# (single-species, dog-only network), per KatarinaTe's request 2026-10-02. Reuses the 6
# already-computed pairs from fig2e_similarity_results_DRAFT.csv unchanged, and only computes
# the 4 new pairs involving dogCD. Same methodology as fig2e_similarity_heatmap_DRAFT.R --
# see that script's header for the full caveats (this is still an unverified reconstruction).

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
existing_results_csv <- "data/results/05_network/pipeline_output/fig2/fig2e_similarity_results_DRAFT.csv"
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

# ---- dogCD term set: same top-community + significance logic as the other 4 groups ----
dogcd_raw <- read_tsv(dogcd_enrichment_file, show_col_types = FALSE)
ccol <- "name.1"
sz <- tapply(suppressWarnings(as.numeric(dogcd_raw$query_size)), dogcd_raw[[ccol]],
             function(v) max(v, na.rm = TRUE))
top_by_size <- names(sz)[which.max(sz)]
nums <- suppressWarnings(as.numeric(str_extract(names(sz), "[0-9]+$")))
stopifnot(identical(top_by_size, names(sz)[which.min(nums)]))  # same sanity check as 01_semantic_clustering.R

dogcd_top <- dogcd_raw %>% filter(.data[[ccol]] == top_by_size) %>%
  transmute(GO_id = native, GO_term = name, p = p_value) %>%
  group_by(GO_id, GO_term) %>% summarize(p = min(p, na.rm = TRUE), .groups = "drop") %>%
  mutate(neglog10p = -log10(p))

groups[["dogCD"]] <- dogcd_top %>% filter(neglog10p >= sig_cutoff) %>% pull(GO_id) %>% unique()
message("dogCD top community: ", top_by_size, " (", nrow(dogcd_top), " terms tested, ",
        length(groups[["dogCD"]]), " significant)")

message("Loading GOSemSim BP semantic data...")
sem_data <- godata(OrgDb = "org.Hs.eg.db", ont = "BP", computeIC = FALSE)
full_bg <- unique(sem_data@geneAnno$GO[sem_data@geneAnno$ONTOLOGY == "BP"])
sim <- function(a, b) mgoSim(a, b, sem_data, measure = "Wang", combine = "BMA")

new_pairs <- tribble(
  ~row,    ~col,
  "dogCD", "OCD",
  "dogCD", "OCD/dogCD",
  "dogCD", "Depression",
  "dogCD", "Schizophrenia"
)

message("Computing 4 new dogCD pairs (", n_perm, " permutations each)...")
new_results <- new_pairs %>%
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
  for (i in seq_len(n_perm)) null[i] <- sim(sample(full_bg, length(a)), sample(full_bg, length(b)))
  null
}

new_results <- new_results %>%
  rowwise() %>%
  mutate(
    null_draws = list(null_mean(a, b)),
    chance = mean(null_draws),
    score = (raw - chance) / (1 - chance),
    pval = (sum(null_draws >= raw) + 1) / (n_perm + 1)
  ) %>%
  ungroup() %>%
  select(row, col, n_overlap, raw, chance, score, pval)

print(new_results)

# ---- combine with the 6 already-computed pairs (unchanged) ----
existing <- read_csv(existing_results_csv, show_col_types = FALSE)
all_results <- bind_rows(existing, new_results)
write_csv(all_results, "fig2e_similarity_results_DRAFT_with_dogCD.csv")

# ---- plot, 5-group upper triangle ----
group_order <- c("dogCD", "OCD", "OCD/dogCD", "Depression", "Schizophrenia")
group_n <- c(dogCD = length(groups[["dogCD"]]), OCD = length(groups[["OCD"]]),
             `OCD/dogCD` = length(groups[["OCD/dogCD"]]), Depression = length(groups[["Depression"]]),
             Schizophrenia = length(groups[["Schizophrenia"]]))
group_labels <- paste0(group_order, "\n(N=", group_n[group_order], ")")
names(group_labels) <- group_order

plot_dat <- all_results %>%
  mutate(
    row = factor(row, levels = rev(group_order[1:4])),
    col = factor(col, levels = group_order[2:5]),
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
  scale_x_discrete(labels = group_labels[group_order[2:5]], position = "bottom") +
  scale_y_discrete(labels = group_labels[rev(group_order[1:4])]) +
  coord_fixed() +
  labs(x = NULL, y = NULL,
       title = "Fig. 2e DRAFT + dogCD -- not verified against source",
       subtitle = "GOSemSim Wang/BMA + 1000-perm chance correction (full GO:BP background)") +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(size = 10, face = "bold", colour = "#8B2E2E"),
    plot.subtitle = element_text(size = 8, colour = "grey30")
  )

ggsave("fig2e_similarity_heatmap_DRAFT_with_dogCD.png", p, width = 6.5, height = 5.5, dpi = 300, bg = "white")
ggsave("fig2e_similarity_heatmap_DRAFT_with_dogCD.pdf", p, width = 6.5, height = 5.5)
message("Wrote fig2e_similarity_heatmap_DRAFT_with_dogCD.png/.pdf and fig2e_similarity_results_DRAFT_with_dogCD.csv")
