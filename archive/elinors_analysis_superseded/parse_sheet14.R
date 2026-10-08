library(tidyverse)

sheet14 <- read_csv(
  "../uploads/gwas_catalog_hits.psych - Sheet14.csv",
  show_col_types = FALSE
)

sheet14_long <- sheet14 %>%
  mutate(gene = str_split(genes, ",\\s*")) %>%
  unnest(gene) %>%
  mutate(gene = str_trim(gene)) %>%
  select(gene, region, P, GWAS)

message(sprintf("Regions: %d -> Genes: %d", nrow(sheet14), nrow(sheet14_long)))
print(sheet14_long)
