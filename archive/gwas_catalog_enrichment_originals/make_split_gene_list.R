library(tidyverse)
library(googlesheets4)

gs4_auth(email = "elinor@broadinstitute.org")

url <- "<Google Sheet link removed>"

SIG_THRESH  <- 4e-7   # significant
SUGG_THRESH <- 1e-6   # suggestive upper bound

# ── 1. Read tabs ─────────────────────────────────────────────────────────────

message("Reading gene list...")
gene_list <- read_sheet(url, sheet = "gene list")
# Preserve original column names; we'll use them as-is for output

message("Reading temp...")
temp_raw <- read_sheet(url, sheet = "temp", col_types = "c")

message("Reading pheno...")
pheno_raw <- read_sheet(url, sheet = "pheno")

# ── 2. Prepare pheno lookup ───────────────────────────────────────────────────

pheno_lookup <- pheno_raw %>%
  rename(code = `Item number`, description = `Short description`) %>%
  filter(!is.na(code), !is.na(description)) %>%
  distinct(code, .keep_all = TRUE) %>%
  select(code, description)

# ── 3. Find the key columns in temp ──────────────────────────────────────────

# The temp tab has two columns named "GWAS" — the 5th column is the phenotype code.
# After reading, googlesheets4 may append "...N" suffixes to duplicates.
message("temp column names: ", paste(names(temp_raw), collapse = " | "))

# Identify by position:
#   col 1: Region_sorted
#   col 3: P
#   col 5: GWAS (phenotype code) — first GWAS column
#   col 19 or similar: genes_251121
p_col     <- names(temp_raw)[3]                           # "P"
gwas_col  <- names(temp_raw)[5]                           # first "GWAS" column
genes_col <- names(temp_raw)[str_detect(names(temp_raw), regex("genes_251121|genes 251121", ignore_case = TRUE))][1]

if (is.na(genes_col)) {
  # fallback: look for "Genes to netcoloc" or similar
  genes_col <- names(temp_raw)[str_detect(names(temp_raw), regex("genes.*netcoloc|netcoloc", ignore_case = TRUE))][1]
}

message("Using: P='", p_col, "', GWAS='", gwas_col, "', genes='", genes_col, "'")

# ── 4. Explode: one row per gene per region-phenotype pair ────────────────────

temp_slim <- temp_raw %>%
  select(P = all_of(p_col), gwas_code = all_of(gwas_col), genes_raw = all_of(genes_col)) %>%
  mutate(P = as.numeric(P)) %>%
  filter(!is.na(P), !is.na(genes_raw),
         genes_raw != "", genes_raw != "no genes", genes_raw != "NA") %>%
  mutate(gene = str_split(genes_raw, ",\\s*")) %>%
  unnest(gene) %>%
  mutate(gene = str_trim(gene)) %>%
  filter(gene != "", gene != "NA") %>%
  left_join(pheno_lookup, by = c("gwas_code" = "code")) %>%
  mutate(pheno_label = coalesce(description, gwas_code))

message("Rows after unnesting: ", nrow(temp_slim), " | Unique genes: ", n_distinct(temp_slim$gene))

# ── 5. Build per-gene phenotype strings ──────────────────────────────────────

# Significant phenotypes: best P per gwas_code ≤ SIG_THRESH
sig_phenos <- temp_slim %>%
  filter(P <= SIG_THRESH) %>%
  group_by(gene, gwas_code, pheno_label) %>%
  summarise(best_P = min(P), .groups = "drop") %>%
  mutate(is_ccdf = str_detect(gwas_code, "^CCDF")) %>%
  arrange(gene, desc(is_ccdf), best_P) %>%
  group_by(gene) %>%
  summarise(sig_pheno_str = paste(unique(pheno_label), collapse = "; "),
            .groups = "drop")

# Identify which gwas_codes are already significant per gene
sig_codes_per_gene <- temp_slim %>%
  filter(P <= SIG_THRESH) %>%
  group_by(gene) %>%
  summarise(sig_codes = list(unique(gwas_code)), .groups = "drop")

# Suggestive-only: > SIG_THRESH AND ≤ SUGG_THRESH AND NOT already significant for this gene
sugg_phenos <- temp_slim %>%
  filter(P > SIG_THRESH, P <= SUGG_THRESH) %>%
  left_join(sig_codes_per_gene, by = "gene") %>%
  mutate(sig_codes = map(sig_codes, ~ if (is.null(.x) || length(.x) == 0) character(0) else as.character(unlist(.x)))) %>%
  filter(!map2_lgl(gwas_code, sig_codes, ~ .x %in% .y)) %>%
  group_by(gene, gwas_code, pheno_label) %>%
  summarise(best_P = min(P), .groups = "drop") %>%
  mutate(is_ccdf = str_detect(gwas_code, "^CCDF")) %>%
  arrange(gene, desc(is_ccdf), best_P) %>%
  group_by(gene) %>%
  summarise(`suggestive dogCD phenotypes` = paste(unique(pheno_label), collapse = "; "),
            .groups = "drop")

message("Genes with significant phenotypes:     ", nrow(sig_phenos))
message("Genes with suggestive-only phenotypes: ", nrow(sugg_phenos))

# ── 6. Join onto gene list ────────────────────────────────────────────────────

new_gene_list <- gene_list %>%
  select(-any_of(c("suggestive dogCD phenotypes", "sig_pheno_str"))) %>%
  left_join(sig_phenos,  by = "gene") %>%
  left_join(sugg_phenos, by = "gene")

# Guarantee column exists even when sugg_phenos has 0 rows
if (!"suggestive dogCD phenotypes" %in% names(new_gene_list)) {
  new_gene_list <- new_gene_list %>% mutate(`suggestive dogCD phenotypes` = NA_character_)
}

new_gene_list <- new_gene_list %>%
  # Overwrite "dogCD phenotypes" with significant-only values
  mutate(`dogCD phenotypes` = sig_pheno_str) %>%
  select(-sig_pheno_str)

# Place "suggestive dogCD phenotypes" immediately after "dogCD phenotypes"
dogcd_col_pos <- which(names(new_gene_list) == "dogCD phenotypes")
if (length(dogcd_col_pos) > 0) {
  before <- names(new_gene_list)[seq_len(dogcd_col_pos)]
  after  <- names(new_gene_list)[seq(dogcd_col_pos + 1, ncol(new_gene_list))]
  after  <- after[after != "suggestive dogCD phenotypes"]
  new_gene_list <- new_gene_list %>%
    select(all_of(before), `suggestive dogCD phenotypes`, all_of(after))
}

message("Output columns: ", paste(names(new_gene_list), collapse = " | "))

# ── 7. Write to Google Sheet ──────────────────────────────────────────────────

message("Writing 'gene list 2'...")
tryCatch(
  write_sheet(new_gene_list, ss = url, sheet = "gene list 2"),
  error = function(e) {
    message("write_sheet failed (", conditionMessage(e), ") — using range_write fallback")
    tryCatch(sheet_add(url, sheet = "gene list 2"), error = function(e2) NULL)
    range_clear(ss = url, sheet = "gene list 2", reformat = FALSE)
    range_write(new_gene_list, ss = url, sheet = "gene list 2",
                range = "A1", col_names = TRUE, reformat = FALSE)
  }
)

message("Done. ", nrow(new_gene_list), " rows written to 'gene list 2'.")
