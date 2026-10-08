#!/usr/bin/env Rscript
# add_venn_region_labels.R -- adds a human-readable "which Fig. 2b Venn region" column per
# comparison to gene_level_network_basis_table.csv (Supplementary Table 8's source data).
#
# Reuses the exact same region definitions plot_venn_hierarchy.R already computes for drawing the
# Venns (region_counts(dep, scz, A), same in_hierarchy_* columns, same 3 comparisons) rather than
# a separate reimplementation, so the new columns are guaranteed consistent with the published
# Fig. 2b panels by construction. Verified below to reproduce the exact same 7 region counts per
# comparison already confirmed against the published figure (see plot_venn_hierarchy.R's header):
#   dogCD      432 517 226 74 16  8  5
#   OCD        399 511 290 71 49 14  8
#   OCD/dogCD  420 511 106 62 28 14 17
#
# Usage: add_venn_region_labels.R <gene_level_network_basis_table.csv> <out.csv>

suppressPackageStartupMessages({
  library(tidyverse)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("Usage: add_venn_region_labels.R <gene_level_network_basis_table.csv> <out.csv>")
gene_table_file <- args[1]
out_file        <- args[2]

# Read twice: `raw` is untouched and is what gets written back out (plus the new columns) --
# pre-existing columns must stay byte-identical, not get silently reformatted (e.g.
# write_csv() would otherwise round-trip "True"/"False" to "TRUE"/"FALSE"). `parsed` is a
# logical-coerced working copy used only to compute region membership.
raw <- read_csv(gene_table_file, col_types = cols(.default = "c")) %>%
  filter(!is.na(gene))
parsed <- raw %>%
  mutate(across(-gene, ~ coalesce(toupper(str_trim(.x)) == "TRUE", FALSE)))

gene_set <- function(col) parsed$gene[parsed[[col]]]

dep <- gene_set("in_hierarchy_Depression")
scz <- gene_set("in_hierarchy_Schizophrenia")

comparisons <- tribble(
  ~column,                  ~label,      ~out_col,
  "in_hierarchy_dogCD",     "dogCD",     "venn_region_dogCD",
  "in_hierarchy_OCD",       "OCD",       "venn_region_OCD",
  "in_hierarchy_OCD/dogCD", "OCD/dogCD", "venn_region_OCD_dogCD"
)

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

# Per-gene region code: which of {Depression, Schizophrenia, <comparison set>} it's in_hierarchy
# for, as a 3-character "DSA" string -- same encoding plot_venn_hierarchy.R uses to pick a shape
# layer/fill colour, here used to assign a label instead. A gene in none of the three (region
# "000") gets NA -- it isn't in any of this particular Venn's three in-hierarchy sets at all.
region_label <- function(region, set_label) {
  case_when(
    region == "100" ~ "Depression only",
    region == "010" ~ "Schizophrenia only",
    region == "001" ~ paste(set_label, "only"),
    region == "110" ~ "Depression + Schizophrenia",
    region == "101" ~ paste("Depression +", set_label),
    region == "011" ~ paste("Schizophrenia +", set_label),
    region == "111" ~ paste("Depression + Schizophrenia +", set_label),
    TRUE ~ NA_character_
  )
}

for (i in seq_len(nrow(comparisons))) {
  column    <- comparisons$column[i]
  set_label <- comparisons$label[i]
  out_col   <- comparisons$out_col[i]
  A <- gene_set(column)

  counts <- region_counts(dep, scz, A)
  stopifnot(sum(counts$n[counts$region %in% c("001", "101", "011", "111")]) == length(A))
  message(set_label, ": ", paste(counts$region, counts$n, sep = "=", collapse = "  "))

  region <- paste0(
    as.integer(parsed$in_hierarchy_Depression), as.integer(parsed$in_hierarchy_Schizophrenia),
    as.integer(parsed[[column]])
  )
  raw[[out_col]] <- region_label(region, set_label)
}

write_csv(raw, out_file)
message("Wrote ", out_file)
