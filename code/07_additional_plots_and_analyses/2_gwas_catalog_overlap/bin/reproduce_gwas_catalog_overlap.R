#!/usr/bin/env Rscript
# Reproduces a collaborator original's parse_gwas_catalog_hits.R join (name/trait matching +
# per-gene psychiatric/OCD/other/temperament hit counts, merged onto the master gene list) without
# needing live Google Sheets access, using local repo-tracked inputs only (moved out of archive/
# into data/ 2026-09-21, since these are real pipeline inputs, not provenance-only material -- see
# REPRODUCIBILITY_AUDIT.md):
#   - data/07_additional_plots_and_analyses/gwas_catalog_hits.csv   (raw GWAS Catalog query results)
#   - data/07_additional_plots_and_analyses/name_conversion_v2.tsv  (local mirror of the "name
#     conversion" sheet that upload_name_conversion.R uploads verbatim)
#   - data/07_additional_plots_and_analyses/gene_list_frozen.csv    (frozen master gene list — see below)
#
# The original scripts' `gene_list` came from a "gene list" Google Sheet tab (region assignment,
# network-community membership, a hand-assigned "rating", free-text "comment") with no file
# equivalent anywhere in this repo and no generating script — it reads as manually curated, not
# computed. `../gene_list_frozen.csv` is a best-effort frozen stand-in for that tab, derived by
# stripping the catalog-derived columns back off the one manual export that WAS already checked
# into this repo (`gene_list_gwas_catalog_overlap_with_region.csv`). Two columns referenced by
# gwas_make_googlesheet.R (`Strom_et_al`, `comment`) are dropped by that same export before it
# ever reaches a file, so they are NOT recoverable from anything in this repo and are simply
# absent from the frozen file. If a real export of the raw "gene list" tab ever becomes available,
# it should replace `gene_list_frozen.csv` directly (same columns, no script changes needed).
#
# Because gene_list_frozen.csv was itself derived from gene_list_gwas_catalog_overlap_with_region.csv,
# reproducing that same file's psychiatric-hit-count columns from it is a consistency check on this
# script's join logic, not independent validation — it confirms the join is wired correctly, not
# that the underlying data is right.
#
# Outputs (written to the working directory):
#   signif_hits_in_gwas_catalog.csv         - one row per gene/study/trait hit that passed name
#                                             conversion (mirrors the "signif hits in gwas catalog" tab)
#   gene_catalog_hit_summary.csv            - one row per gene: disorder/other/temperament hit
#                                             counts only (no gene_list columns)
#   gene_list_with_catalog_hit_counts.csv   - gene_list_frozen.csv's columns with the above hit
#                                             counts left-joined on (mirrors the mechanical part of
#                                             "gene list w gwas catalog overlap raw"; does not
#                                             attempt the frozen CSV's hand-formatted
#                                             "Depression or schizophrenia?"/pubmed-string columns,
#                                             which involve formatting choices this script can't
#                                             verify without live execution)
#   gene_list_overlap_crosscheck.csv        - the frozen CSV's own psychiatric hit count next to
#                                             this script's recomputed count, for genes in both,
#                                             flagging any mismatch
#
# Takes 4 positional args: catalog_hits_file, name_conversion_file, frozen_overlap_file,
# frozen_gene_list_file (paths staged by Nextflow; see process_definitions.nf).

suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
catalog_hits_file <- args[1]
name_conversion_file <- args[2]
frozen_overlap_file <- args[3]
frozen_gene_list_file <- args[4]

df_raw <- as_tibble(read.csv(catalog_hits_file, header = TRUE)) %>%
  mutate(pvalue_computed = pvalue_mantissa * 10^pvalue_exponent)

convert <- read_tsv(name_conversion_file, show_col_types = FALSE) %>%
  mutate(include = as.logical(include)) %>%
  filter(include) %>%
  mutate(
    disorder_label = if_else(
      short %in% c("Depression", "Schizophrenia", "Obsessive-compulsive disorder"),
      short, ""
    )
  )

missing_traits <- setdiff(unique(df_raw$reported_trait), convert$reported_trait)
if (length(missing_traits) > 0) {
  warning(
    length(missing_traits), " reported_trait(s) in gwas_catalog_hits.csv are not in ",
    "name_conversion_v2.tsv and are dropped by the inner join below:\n",
    paste(sort(missing_traits), collapse = "\n")
  )
}

df_hits <- df_raw %>%
  select(-any_of(c("disorder", "is_psychiatric", "is_temperament", "temperament_category"))) %>%
  inner_join(convert %>% select(reported_trait, short, category, disorder_label), by = "reported_trait") %>%
  arrange(queried_gene, pvalue_computed)

write_csv(df_hits, "signif_hits_in_gwas_catalog.csv")

summarize_category <- function(df_cat, count_col_name) {
  df_cat %>%
    group_by(queried_gene, short) %>%
    summarize(pvalue_computed = min(pvalue_computed), n_assoc = n(), .groups = "drop") %>%
    group_by(queried_gene) %>%
    summarize("{count_col_name}" := sum(n_assoc), .groups = "drop")
}

disorder_counts <- summarize_category(df_hits %>% filter(category == "Disorder"), "n_disorder_assoc")
other_counts <- summarize_category(df_hits %>% filter(category == "Other"), "n_other_assoc")
temperament_counts <- summarize_category(df_hits %>% filter(category == "Temperament"), "n_temperament_assoc")

gene_summary <- disorder_counts %>%
  full_join(other_counts, by = "queried_gene") %>%
  full_join(temperament_counts, by = "queried_gene") %>%
  mutate(across(starts_with("n_"), ~ replace_na(., 0))) %>%
  arrange(queried_gene)

write_csv(gene_summary, "gene_catalog_hit_summary.csv")

if (file.exists(frozen_gene_list_file)) {
  gene_list <- read_csv(frozen_gene_list_file, show_col_types = FALSE)

  gene_list_with_hits <- gene_list %>%
    left_join(gene_summary, by = c("gene" = "queried_gene")) %>%
    mutate(across(starts_with("n_"), ~ replace_na(., 0)))

  write_csv(gene_list_with_hits, "gene_list_with_catalog_hit_counts.csv")
  message("Wrote gene_list_with_catalog_hit_counts.csv (", nrow(gene_list_with_hits), " genes).")
} else {
  message("Frozen gene list not found at ", frozen_gene_list_file, "; skipping that join.")
}

if (file.exists(frozen_overlap_file)) {
  frozen <- read_csv(frozen_overlap_file, show_col_types = FALSE) %>%
    select(gene, frozen_n_psych = `num pysch assoc in catalog`) %>%
    filter(!is.na(gene))

  crosscheck <- frozen %>%
    inner_join(
      gene_summary %>% mutate(recomputed_n_psych = n_disorder_assoc + n_temperament_assoc) %>%
        select(queried_gene, recomputed_n_psych),
      by = c("gene" = "queried_gene")
    ) %>%
    mutate(matches = coalesce(frozen_n_psych, 0) == recomputed_n_psych)

  write_csv(crosscheck, "gene_list_overlap_crosscheck.csv")

  n_mismatch <- sum(!crosscheck$matches)
  message(
    "Cross-check vs. frozen gene_list_gwas_catalog_overlap_with_region.csv: ",
    nrow(crosscheck), " genes compared, ", n_mismatch, " mismatch(es)."
  )
  if (n_mismatch > 0) {
    message(
      "Mismatches are expected wherever the frozen file's psychiatric count also depends on the ",
      "manually curated gene_list sheet's own row selection — see this script's header."
    )
  }
} else {
  message("Frozen overlap file not found at ", frozen_overlap_file, "; skipping cross-check.")
}

message("Done. Wrote signif_hits_in_gwas_catalog.csv (", nrow(df_hits), " rows) and ",
        "gene_catalog_hit_summary.csv (", nrow(gene_summary), " genes).")
