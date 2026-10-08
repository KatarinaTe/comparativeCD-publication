#!/usr/bin/env Rscript
library(tidyverse)
library(scales)
library(stringr)
library(tidytext)
library(cowplot)

#========================
# Settings
#========================

plot_ds          <- "OCD_CCD"
communities      <- c("C185", "C186", "C197","C187","C193")

pvalue_max       <- 1e-5
max_sets         <- 10
min_intersection <- 5
min_term_size    <- 5

plot_width       <- 7
plot_height      <- 5
plot_dpi         <- 300

cache_file       <- "all_go_combined.rds"

point_alpha      <- 0.95
segment_alpha    <- 0.8
segment_width    <- 0.4
gene_label_size  <- 2.3

color_values <- c("#2C7BB6", "#7B3294", "#D7191C")

input_files <- tibble(
  ds = c("DEP_only", "DEP_CCD", "OCD_only", "OCD_CCD", "SCH_only", "SCH_CCD"),
  file = c(
    "../go_enrichment_tsvs/DEP_ONLY_hierarchy_full_GO_enrichment.tsv",
    "../go_enrichment_tsvs/DEP_CCD_hierachy_full_GO_enrichment.tsv",
    "../go_enrichment_tsvs/OCD_ONLY_hierachy_full_GO_enrichment.tsv",
    "../go_enrichment_tsvs/Compulsive_hierachy_full_GO_enrichment_251023.tsv",
    "../go_enrichment_tsvs/SCH_ONLY_hierachy_full_GO_enrichment.tsv",
    "../go_enrichment_tsvs/SCH_CCD_hierachy_full_GO_enrichment.tsv"
  )
)

#========================
# Helpers
#========================

wrap_one_newline_middle <- function(x) {
  max_width <- max(nchar(x), na.rm = TRUE)
  cutoff <- ceiling(max_width / 2)
  
  purrr::map_chr(x, function(s) {
    n <- nchar(s)
    
    if (is.na(s) || n <= cutoff) {
      return(s)
    }
    
    spaces <- gregexpr(" ", s, fixed = TRUE)[[1]]
    
    if (length(spaces) == 1 && spaces[1] == -1) {
      return(s)
    }
    
    mid <- n / 2
    split_at <- spaces[which.min(abs(spaces - mid))]
    
    paste0(
      substr(s, 1, split_at - 1),
      "\n",
      substr(s, split_at + 1, n)
    )
  })
}

make_outfile <- function(ext) {
  paste0(
    "GO_enrichment_",
    plot_ds,
    "_facets.p",
    -log10(pvalue_max),
    ".top",
    max_sets,
    ".",
    ext
  )
}

#========================
# Read / cache concatenated input files
#========================

cache_is_stale <-
  !file.exists(cache_file) ||
  any(file.info(input_files$file)$mtime > file.info(cache_file)$mtime)

if (!cache_is_stale) {
  
  message("Loading cached GO table: ", cache_file)
  all_go_raw <- readRDS(cache_file)
  
} else {
  
  message("Reading GO enrichment files")
  
  all_go_raw <- input_files %>%
    mutate(
      file_id = basename(file) %>% str_remove("\\.tsv$"),
      dat = purrr::map(file, ~ readr::read_tsv(.x, show_col_types = FALSE))
    ) %>%
    tidyr::unnest(dat) %>% select(-...1)
  
  all_go_raw <- all_go_raw %>%
    mutate(
      community = if_else(
        is.na(community) & !is.na(name.1),
        name.1,
        community)) %>% 
      select(-name.1)
  
  saveRDS(
    all_go_raw,
    cache_file,
    compress = "xz"
  )
  
  message("Cached GO table: ", cache_file)
}

#========================
# Clean / calculate enrichment values
#========================

all_go <- all_go_raw %>%
  mutate(
    expected = (query_size * term_size) / effective_domain_size,
    fold_enrichment = intersection_size / expected
  ) %>%
  filter(term_size > min_term_size)

#========================
# Plot function
#========================

plot_go_enrichment_facets <- function(all_go, in_ds, communities) {
  
  pd <- all_go %>%
    filter(
      ds == in_ds,
      community %in% communities,
      p_value < pvalue_max,
      intersection_size > min_intersection
    ) %>%
    mutate(
      term_size = round(term_size),
      neg_log10_p = -log10(p_value),
      y_lab_raw = paste0(native, ": ", name, " (", term_size, ")")
    ) %>%
    group_by(community) %>%
    arrange(desc(fold_enrichment), .by_group = TRUE) %>%
    slice_head(n = max_sets) %>%
    ungroup() %>%
    mutate(
      y_lab = wrap_one_newline_middle(y_lab_raw),
      y_lab = tidytext::reorder_within(y_lab, fold_enrichment, community),
      gene_label = as.character(intersection_size),
      label_x = fold_enrichment + 0.04 * max(fold_enrichment, na.rm = TRUE)
    )
  
  message(paste(pd %>% dplyr::select(community,native) %>% distinct(),sep="\n"))
  p <- ggplot(
    pd,
    aes(
      x = fold_enrichment,
      y = y_lab,
      colour = neg_log10_p,
      size = intersection_size
    )
  )
  
  p <- p + geom_segment(
    aes(x = 0, xend = fold_enrichment, yend = y_lab),
    linewidth = segment_width,
    alpha = segment_alpha
  )
  
  p <- p + geom_point(alpha = point_alpha)
  
  p <- p + geom_text(
    aes(x = label_x, label = gene_label),
    size = gene_label_size,
    colour = "grey20",
    hjust = 0,
    show.legend = FALSE
  )
  
  p <- p + facet_grid(
    community ~ .,
    scales = "free_y",
    space = "free_y"
  )
  
  p <- p + tidytext::scale_y_reordered()
  
  p <- p + scale_x_continuous(
    expand = expansion(mult = c(0.02, 0.18))
  )
  
  p <- p + scale_color_gradientn(
    colours = color_values,
    name = expression(-log[10](p))
  )
  
  p <- p + scale_size_continuous(
    name = "Genes",
    range = c(1.5, 5)
  )
  
  p <- p + labs(
    x = "Fold enrichment",
    y = "GO ID: GO term (term size)"
  )
  
  p <- p + theme_cowplot()
  
  p <- p + theme(
    strip.background = element_rect(
      fill = "grey95",
      colour = "grey70",
      linewidth = 0.3
    ),
    strip.text.y = element_text(
      size = 9,
      face = "bold",
      angle = 0,
      margin = margin(2, 2, 2, 2)
    ),
    panel.spacing.y = unit(0.8, "lines"),
    panel.border = element_rect(
      colour = "grey80",
      fill = NA,
      linewidth = 0.3
    ),
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(size = 8),
    axis.title.y = element_text(size = 9),
    axis.title.x = element_text(size = 9),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7)
  )
  
  p
}

#========================
# Run
#========================

p <- plot_go_enrichment_facets(
  all_go = all_go,
  in_ds = plot_ds,
  communities = communities
)

p

ggsave(
  filename = make_outfile("pdf"),
  plot = p,
  width = plot_width,
  height = plot_height,
  units = "in"
)

ggsave(
  filename = make_outfile("png"),
  plot = p,
  width = plot_width,
  height = plot_height,
  units = "in",
  dpi = plot_dpi
)
