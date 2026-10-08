library(tidyverse)
library(googlesheets4)

# Helper to remove underscores from column names
clean_names <- function(df) {
  names(df) <- str_replace_all(names(df), "_", " ")
  names(df) <- str_replace_all(names(df), "psychiatric", "psych")
  names(df) <- str_replace_all(names(df), "nonpsych", "non psych") 
  names(df) <- str_replace_all(names(df), "psych score","% hits for psych vs non psych")
  df
}
# Read in data
results <- as_tibble(read_csv("gwas_catalog_hits.csv")) %>%
  mutate(pubmed_id = as.character(pubmed_id)) %>% 
  mutate(queried_gene=if_else(queried_gene=="NALF1","FAM155A",queried_gene))


# Read in updated psychiatric labels
updated_labels <- read_csv("updated_is_psychiatric_labels.csv") %>%
  select(reported_trait, pubmed_id, is_neuropsych) %>%
  mutate(pubmed_id = as.character(pubmed_id))

# Join and update is_psychiatric
results <- results %>%
  left_join(updated_labels, by = c("reported_trait", "pubmed_id")) %>%
  mutate(is_psychiatric = !is.na(is_neuropsych) & is_neuropsych) %>%
  select(-is_neuropsych)






# Calculate gene scores
# calculate gene scores as before but without psychiatric_disorders column
gene_scores <- results %>%
  group_by(queried_gene, reported_trait, is_psychiatric, pubmed_id, disorder) %>%
  summarise(
    min_p = min(pvalue),
    .groups = "drop"
  ) %>%
  mutate(is_significant = min_p < 5e-8) %>%
  group_by(queried_gene) %>%
  summarise(
    n_psychiatric_traits    = sum(is_significant & is_psychiatric),
    n_nonpsychiatric_traits = sum(is_significant & !is_psychiatric),
    psychiatric_score       = n_psychiatric_traits /
      (n_psychiatric_traits + n_nonpsychiatric_traits),
    # top psychiatric hit
    top_psych_pvalue        = if (any(is_psychiatric)) min(min_p[is_psychiatric]) else NA_real_,
    top_psych_trait         = if (any(is_psychiatric)) reported_trait[is_psychiatric][which.min(min_p[is_psychiatric])] else NA_character_,
    top_psych_pmid          = if (any(is_psychiatric)) pubmed_id[is_psychiatric][which.min(min_p[is_psychiatric])] else NA_character_,
    .groups = "drop"
  )

sheet_url <- "<Google Sheet link removed>"

# Read in gene list sheet
gene_list <- read_sheet(sheet_url, sheet = "gene list") %>%
  rename_with(~ str_replace_all(., " ", "_"))

# --- Sheet 1: gene scores ---
merged_scores <- gene_list %>% select(-dogCD_P,-Strom_et_al,-comment) %>% 
  left_join(gene_scores, by = c("gene" = "queried_gene")) %>%
  arrange(-psychiatric_score) 

write_sheet(
  data  = clean_names(merged_scores),
  ss    = sheet_url,
  sheet = "scores by gene"
)
message("Sheet 1 done — ", nrow(merged_scores), " genes written to 'gene scores'.")

# --- Sheet 2: psychiatric catalog hits ---
catalog <- results %>%
  filter(is_psychiatric) %>%
  select(
    gene        = queried_gene,
    pvalue,
    reported_trait,
    pubmed_id,
    pub_date,
    title,
    disorder,
    sample_size
  ) %>%
  distinct()

merged_catalog <- gene_list %>% select(gene,dogCD_P,rating) %>% 
  left_join(catalog, by = "gene") %>% select(-disorder) %>% distinct() %>% 
  arrange(pvalue)

write_sheet(
  data  = clean_names(merged_catalog),
  ss    = sheet_url,
  sheet = "psychiatric hits"
)
message("Sheet 2 done — ", nrow(merged_catalog), " rows written to 'psychiatric hits'.")

# --- Sheet 3: top 5 psychiatric traits per gene ---
top5 <- gene_list %>% select(-Strom_et_al,-comment) %>% 
  inner_join(catalog, by = "gene") %>%
  group_by(gene) %>%
  slice_min(order_by = pvalue, n = 5, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(desc(dog_seed),rating,dogCD_P,gene, pvalue)

write_sheet(
  data  = clean_names(top5),
  ss    = sheet_url,
  sheet = "top 5 traits"
)
message("Sheet 3 done — ", nrow(top5), " rows written to 'top 5 traits'.")