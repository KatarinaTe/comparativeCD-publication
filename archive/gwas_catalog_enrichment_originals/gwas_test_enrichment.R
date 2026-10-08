# ============================================================
# Gene enrichment for psychiatric / OCD GWAS Catalog traits
# across multiple GWAS p-value thresholds
#
# Input: gwas_catalog_hits.csv
# Outputs:
#   gwas_gene_psych_enrichment_by_threshold.csv
#   gwas_gene_trait_enrichment_summary_by_threshold.csv
# ============================================================

library(dplyr)
library(readr)
library(tibble)
library(stringr)
library(purrr)

in_file <- "gwas_catalog_hits.csv"

gwas_thresholds <- c(5e-8, 1e-10, 1e-12, 1e-20)

out_psych <- "gwas_gene_psych_enrichment_by_threshold.csv"
out_summary <- "gwas_gene_trait_enrichment_summary_by_threshold.csv"

dat <- read_csv(in_file, show_col_types = FALSE)

required_cols <- c(
  "queried_gene",
  "study_id",
  "reported_trait",
  "pvalue",
  "is_psychiatric",
  "disorder"
)

missing_cols <- setdiff(required_cols, names(dat))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

dat <- dat %>%
  mutate(
    queried_gene = as.character(queried_gene),
    study_id = as.character(study_id),
    reported_trait = as.character(reported_trait),
    pvalue = as.numeric(pvalue),
    is_psychiatric = case_when(
      is_psychiatric %in% TRUE ~ TRUE,
      is_psychiatric %in% FALSE ~ FALSE,
      str_to_lower(as.character(is_psychiatric)) == "true" ~ TRUE,
      str_to_lower(as.character(is_psychiatric)) == "false" ~ FALSE,
      TRUE ~ NA
    ),
    disorder = as.character(disorder),
    disorder = na_if(disorder, "NA"),
    disorder = na_if(disorder, "NULL"),
    disorder = na_if(disorder, "null"),
    disorder = na_if(disorder, "")
  )

run_gene_enrichment <- function(dat_sub, flag_col, label, threshold) {
  
  genes <- sort(unique(dat_sub$queried_gene))
  
  map_dfr(genes, function(gene) {
    
    in_gene <- dat_sub$queried_gene == gene
    is_hit <- dat_sub[[flag_col]] %in% TRUE
    
    a <- sum(in_gene & is_hit, na.rm = TRUE)
    b <- sum(in_gene & !is_hit, na.rm = TRUE)
    c <- sum(!in_gene & is_hit, na.rm = TRUE)
    d <- sum(!in_gene & !is_hit, na.rm = TRUE)
    
    mat <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
    
    ft_greater <- fisher.test(mat, alternative = "greater")
    ft_two_sided <- fisher.test(mat, alternative = "two.sided")
    
    tibble(
      gwas_p_threshold = threshold,
      queried_gene = gene,
      test = label,
      n_gene_hits = a + b,
      n_gene_target_hits = a,
      n_gene_other_hits = b,
      n_background_target_hits = c,
      n_background_other_hits = d,
      prop_gene_target = ifelse((a + b) > 0, a / (a + b), NA_real_),
      prop_background_target = ifelse((c + d) > 0, c / (c + d), NA_real_),
      odds_ratio = unname(ft_greater$estimate),
      p_enrichment = ft_greater$p.value,
      p_two_sided = ft_two_sided$p.value
    )
  }) %>%
    mutate(
      p_enrichment_fdr = p.adjust(p_enrichment, method = "BH"),
      p_two_sided_fdr = p.adjust(p_two_sided, method = "BH")
    ) %>%
    arrange(gwas_p_threshold, p_enrichment, desc(odds_ratio))
}

run_threshold <- function(threshold) {
  
  message("\nRunning threshold: p <= ", threshold)
  
  dat_sub <- dat %>%
    filter(
      !is.na(pvalue),
      pvalue <= threshold,
      !is.na(queried_gene),
      queried_gene != "",
      !is.na(study_id),
      study_id != "",
      !is.na(reported_trait),
      reported_trait != "",
      !is.na(is_psychiatric)
    ) %>%
    distinct(queried_gene, study_id, reported_trait, .keep_all = TRUE) %>%
    mutate(
      psych_flag = is_psychiatric,
      ocd_flag = disorder == "OCD"
    )
  
  message("  Deduplicated study-trait hits: ", nrow(dat_sub))
  message("  Genes tested: ", n_distinct(dat_sub$queried_gene))
  message("  Psychiatric hits: ", sum(dat_sub$psych_flag, na.rm = TRUE))
  message("  OCD hits: ", sum(dat_sub$ocd_flag, na.rm = TRUE))
  
  if (nrow(dat_sub) == 0 || n_distinct(dat_sub$queried_gene) < 2) {
    warning("Skipping threshold ", threshold, ": not enough data.")
    return(NULL)
  }
  
  psych <- run_gene_enrichment(dat_sub, "psych_flag", "psychiatric", threshold)
  ocd <- run_gene_enrichment(dat_sub, "ocd_flag", "OCD", threshold)
  
  summary <- dat_sub %>%
    group_by(queried_gene) %>%
    summarize(
      gwas_p_threshold = threshold,
      n_study_trait_hits = n(),
      n_psych_hits = sum(psych_flag, na.rm = TRUE),
      n_ocd_hits = sum(ocd_flag, na.rm = TRUE),
      n_depression_hits = sum(disorder == "depression", na.rm = TRUE),
      n_schizophrenia_hits = sum(disorder == "schizophrenia", na.rm = TRUE),
      n_other_psych_hits = sum(disorder == "other_psychiatric", na.rm = TRUE),
      prop_psych = n_psych_hits / n_study_trait_hits,
      prop_ocd = n_ocd_hits / n_study_trait_hits,
      .groups = "drop"
    ) %>%
    left_join(
      psych %>%
        select(
          gwas_p_threshold,
          queried_gene,
          psych_odds_ratio = odds_ratio,
          psych_p_enrichment = p_enrichment,
          psych_p_enrichment_fdr = p_enrichment_fdr
        ),
      by = c("gwas_p_threshold", "queried_gene")
    ) %>%
    left_join(
      ocd %>%
        select(
          gwas_p_threshold,
          queried_gene,
          ocd_odds_ratio = odds_ratio,
          ocd_p_enrichment = p_enrichment,
          ocd_p_enrichment_fdr = p_enrichment_fdr
        ),
      by = c("gwas_p_threshold", "queried_gene")
    )
  
  list(psych = psych, ocd = ocd, summary = summary)
}

threshold_results <- map(gwas_thresholds, run_threshold)
threshold_results <- compact(threshold_results)

psych_enrichment <- map_dfr(threshold_results, "psych")
gene_summary <- map_dfr(threshold_results, "summary")

write_csv(psych_enrichment, out_psych)
write_csv(gene_summary, out_summary)

message("\nWritten: ", out_psych)
message("Written: ", out_ocd)
message("Written: ", out_summary)

message("\nTop psychiatric-enriched genes by threshold:")
print(
  psych_enrichment %>%
    group_by(gwas_p_threshold) %>%
    slice_min(p_enrichment, n = 10, with_ties = FALSE) %>%
    ungroup() %>%
    select(
      gwas_p_threshold,
      queried_gene,
      n_gene_hits,
      n_gene_target_hits,
      prop_gene_target,
      prop_background_target,
      odds_ratio,
      p_enrichment,
      p_enrichment_fdr
    )
)
