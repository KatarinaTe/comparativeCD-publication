library(tidyverse)
library(googlesheets4)

gs4_auth(email = "elinor@broadinstitute.org")

# safe_write_sheet() fails with "Properties not recognized: tables" on Sheets that
# use Google's newer Tables feature.  This wrapper catches that specific error
# and falls back to range_clear + range_write, which skips the schema check.
safe_write_sheet <- function(data, ss, sheet) {
  tryCatch(
    write_sheet(data, ss = ss, sheet = sheet),
    error = function(e) {
      if (!grepl("tables", conditionMessage(e), ignore.case = TRUE)) stop(e)
      message("  [safe_write_sheet] 'tables' schema error on sheet '", sheet,
              "' — falling back to range_write")
      range_clear(ss = ss, sheet = sheet, reformat = FALSE)
      range_write(data, ss = ss, sheet = sheet,
                  range = "A1", col_names = TRUE, reformat = FALSE)
    }
  )
}

#========================
# Settings
#========================

url    <- "<Google Sheet link removed>"
infile <- "gwas_catalog_hits.csv"

#========================
# Load and compute p-value from exponent
#========================

df_raw <- as_tibble(read.csv(infile, header = TRUE)) %>%
  mutate(
    pvalue_computed = pvalue_mantissa * 10^pvalue_exponent
  )

#========================
# Name conversion table (maps reported_trait -> short name + category)
#========================

convert <- as_tibble(read_sheet(url, "name conversion")) %>%
  filter(include) %>%
  mutate(
    disorder_label = if_else(
      short %in% c("Depression", "Schizophrenia", "Obsessive-compulsive disorder"),
      short, ""
    )
  )

#========================
# Check for traits in df_raw not covered by name conversion
#========================

all_convert_traits <- as_tibble(read_sheet(url, "name conversion")) %>% pull(reported_trait)
missing_traits <- setdiff(unique(df_raw$reported_trait), all_convert_traits)

if (length(missing_traits) > 0) {
  stop(
    length(missing_traits), " reported_trait(s) in ", infile,
    " are not in 'name conversion' and will be silently dropped.\n",
    "Add them to the sheet before re-running.\n\n",
    "Missing traits:\n", paste(sort(missing_traits), collapse = "\n")
  )
}

#========================
# All significant hits passing name conversion filter
#========================

df_hits <- df_raw %>%
  select(-any_of(c("disorder", "is_psychiatric", "is_temperament", "temperament_category"))) %>%
  inner_join(convert %>% select(reported_trait, short, category, disorder_label),
             by = "reported_trait") %>%
  arrange(queried_gene, pvalue_computed)

safe_write_sheet(df_hits, ss = url, sheet = "signif hits in gwas catalog")

#========================
# Summary: Disorder hits per gene
#========================

df_disorder <- df_hits %>% filter(category == "Disorder")

small_psych <- df_disorder %>%
  group_by(queried_gene, short, category) %>%
  summarize(pvalue_computed = min(pvalue_computed), n_assoc = n(), .groups = "drop") %>%
  left_join(df_disorder %>% select(queried_gene, short, pubmed_id) %>% distinct(),
            by = c("queried_gene", "short")) %>%
  arrange(pubmed_id) %>%
  group_by(queried_gene, short, category, pvalue_computed, n_assoc) %>%
  summarize(pubmed_id = paste(pubmed_id, collapse = "; "), n_pmids = n(),
            .groups = "drop") %>%
  mutate(
    dep = if_else(short == "Depression", 1L, 0L),
    sch = if_else(short == "Schizophrenia", 1L, 0L),
    pubmed_id = if_else(
      n_pmids == 1,
      paste0(short, " (PMID ", pubmed_id, ")"),
      paste0(short, " (PMIDs ", pubmed_id, ")")
    )
  ) %>%
  arrange(desc(dep + sch), pvalue_computed) %>%
  group_by(queried_gene) %>%
  summarize(
    catalog_disorder      = if_else(
      max(dep) > 0,
      if_else(max(sch) > 0, "Both", "Depression"),
      if_else(max(sch) > 0, "Schizophrenia", "")
    ),
    n_disorder_assoc      = sum(n_assoc),
    gwas_catalog_disorder = paste(pubmed_id, collapse = "; "),
    .groups = "drop"
  )

#========================
# Summary: Other (neuroimaging, cognitive, sleep, etc.) hits per gene
#========================

df_other <- df_hits %>% filter(category == "Other")

small_other <- df_other %>%
  group_by(queried_gene, short) %>%
  summarize(pvalue_computed = min(pvalue_computed), n_assoc = n(), .groups = "drop") %>%
  left_join(df_other %>% select(queried_gene, short, pubmed_id) %>% distinct(),
            by = c("queried_gene", "short")) %>%
  arrange(pubmed_id) %>%
  group_by(queried_gene, short, pvalue_computed, n_assoc) %>%
  summarize(pubmed_id = paste(pubmed_id, collapse = "; "), n_pmids = n(),
            .groups = "drop") %>%
  mutate(
    pubmed_id = if_else(
      n_pmids == 1,
      paste0(short, " (PMID ", pubmed_id, ")"),
      paste0(short, " (PMIDs ", pubmed_id, ")")
    )
  ) %>%
  arrange(pvalue_computed) %>%
  group_by(queried_gene) %>%
  summarize(
    n_other_assoc      = sum(n_assoc),
    gwas_catalog_other = paste(pubmed_id, collapse = "; "),
    .groups = "drop"
  )

#========================
# Summary: Temperament hits per gene
#========================

df_temp <- df_hits %>% filter(category == "Temperament")

small_temp <- df_temp %>%
  group_by(queried_gene, short) %>%
  summarize(pvalue_computed = min(pvalue_computed), n_assoc = n(), .groups = "drop") %>%
  left_join(df_temp %>% select(queried_gene, short, pubmed_id) %>% distinct(),
            by = c("queried_gene", "short")) %>%
  arrange(pubmed_id) %>%
  group_by(queried_gene, short, pvalue_computed, n_assoc) %>%
  summarize(pubmed_id = paste(pubmed_id, collapse = "; "), n_pmids = n(),
            .groups = "drop") %>%
  mutate(
    pubmed_id = if_else(
      n_pmids == 1,
      paste0(short, " (PMID ", pubmed_id, ")"),
      paste0(short, " (PMIDs ", pubmed_id, ")")
    )
  ) %>%
  arrange(pvalue_computed) %>%
  group_by(queried_gene) %>%
  summarize(
    n_temperament_assoc      = sum(n_assoc),
    temperament_categories   = paste(sort(unique(short)), collapse = "; "),
    gwas_catalog_temperament = paste(pubmed_id, collapse = "; "),
    .groups = "drop"
  )

#========================
# Merge all summaries onto gene list
#========================

gene_list <- read_sheet(url, sheet = "gene list")

small <- gene_list %>%
  full_join(small_psych %>% rename(gene = queried_gene), by = "gene") %>%
  full_join(small_other %>% rename(gene = queried_gene), by = "gene") %>%
  full_join(small_temp  %>% rename(gene = queried_gene), by = "gene") %>%
  arrange(gene) %>%
  rename(
    `dog GWAS gene`                = dog_gene,
    `dog seed`                     = dog_seed,
    `human seed`                   = human_seed,
    `OCD dogCD systems map`        = OCD_dogCD,
    `shared across all maps`       = shared_across_all,
    `dogCD P`                      = `dogCD P`,
    `n psych disorder assoc`        = n_disorder_assoc,
    `gwas catalog (psych disorder)` = gwas_catalog_disorder,
    `n other assoc`                = n_other_assoc,
    `gwas catalog (other)`         = gwas_catalog_other,
    `n temperament assoc`          = n_temperament_assoc,
    `gwas catalog (temperament)`   = gwas_catalog_temperament
  ) %>%
  select(-any_of(c("catalog_disorder", "temperament_categories"))) %>%
  mutate(across(c(`n psych disorder assoc`, `n other assoc`, `n temperament assoc`),
                ~ replace_na(., 0)),
         `n dogCD phenotypes` = if_else(
           is.na(`dogCD phenotypes`) | `dogCD phenotypes` == "",
           0L,
           str_count(`dogCD phenotypes`, ";") + 1L
         )) %>%
  relocate(`n dogCD phenotypes`, .before = `dogCD phenotypes`) %>%
  relocate(`n psych disorder assoc`, `n temperament assoc`, `n other assoc`,
           `gwas catalog (psych disorder)`, `gwas catalog (temperament)`, `gwas catalog (other)`,
           .after = last_col())

names(small) <- str_replace_all(names(small), "_", " ")

safe_write_sheet(small, ss = url, sheet = "gene list w gwas catalog overlap raw")

message("Done. Total hits: ", nrow(df_hits),
        " | Disorder: ", nrow(df_disorder),
        " | Other: ", nrow(df_other),
        " | Temperament: ", nrow(df_temp))
