library(googlesheets4)
library(dplyr)
library(tidyr)

SHEET_ID <- "<Google Sheet link removed>"

# Info columns per sheet: gene metadata paired with the first gene-set column
INFO_COLS <- list(
  "Supplementary Table 8" = c("Seed genes", "Shared/not-shared", "Illustration"),
  "Supplementary Table 9" = character(0)
)

pivot_gene_sets <- function(sheet_id, sheet_name) {
  raw        <- read_sheet(sheet_id, sheet = sheet_name, col_names = TRUE)
  info_cols  <- INFO_COLS[[sheet_name]]
  anchor_col <- names(raw)[1]   # first column is always the gene-info anchor

  # ── Extract gene metadata (anchor gene → info columns) ───────────────────
  if (length(info_cols) > 0) {
    gene_info <- raw %>%
      select(gene = all_of(anchor_col), all_of(info_cols)) %>%
      filter(!is.na(gene), !gene %in% c("#N/A", "")) %>%
      mutate(gene = trimws(as.character(gene)))
  }

  # ── Pivot gene-set columns wide (drop info cols first) ───────────────────
  wide <- raw %>%
    select(-any_of(info_cols)) %>%
    pivot_longer(everything(), names_to = "set", values_to = "gene") %>%
    filter(!is.na(gene), !gene %in% c("#N/A", "")) %>%
    mutate(gene = trimws(as.character(gene)), present = TRUE) %>%
    pivot_wider(names_from = set, values_from = present, values_fill = FALSE) %>%
    arrange(gene)

  # ── Merge gene info back ─────────────────────────────────────────────────
  if (length(info_cols) > 0) {
    wide <- wide %>% left_join(gene_info, by = "gene")
  }

  wide
}

write_sheet_safe <- function(data, ss, sheet) {
  if (sheet %in% sheet_names(ss)) sheet_write(data, ss, sheet)
  else { sheet_add(ss, sheet = sheet); sheet_write(data, ss, sheet) }
  cat(sprintf("Written '%s': %d genes x %d cols\n", sheet, nrow(data), ncol(data) - 1))
}

wide8 <- pivot_gene_sets(SHEET_ID, "Supplementary Table 8")
wide9 <- pivot_gene_sets(SHEET_ID, "Supplementary Table 9")

# ── Merge "extra" into ST8 ───────────────────────────────────────────────────
extra <- read_sheet(SHEET_ID, sheet = "extra", col_names = TRUE)
gene_col_extra <- names(extra)[1]   # assume first column is the gene key
extra <- extra %>% rename(gene = all_of(gene_col_extra))

# Left-join extra info onto ST8 genes only
wide8_merged <- wide8 %>%
  left_join(
    extra %>% select(gene, C185, C186, C197),
    by = "gene"
  ) %>%
  arrange(gene)

write_sheet_safe(wide8,        SHEET_ID, "ST8 wide")
write_sheet_safe(wide8_merged, SHEET_ID, "ST8 wide + extra")
write_sheet_safe(wide9,        SHEET_ID, "ST9 wide")
