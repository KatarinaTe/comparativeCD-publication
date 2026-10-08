## ============================================================================
## enrichment.R
##
## (1) Build a gene x category table for synaptic genes:
##       - glutamatergic-synapse flag from GO (GO:0098978 + offspring)
##       - pre-/post-synaptic location + functional bucket from SynGO
##       Buckets: presynaptic_release, postsynaptic_receptor, scaffold,
##                signalling_mediator
## (2) Read query gene sets from a Google Sheet: one gene column + any number
##     of TRUE/FALSE columns. Each TRUE/FALSE column defines a query set
##     (the genes flagged TRUE). Run over-representation analysis (Fisher,
##     one-sided "greater") of every query set against every bucket category.
##
## ---------------------------------------------------------------------------
## METHODOLOGY NOTES (matter for a methods section / reviewers)
##  * Background choice DRIVES the result. Modes (BACKGROUND_MODE):
##      "sheet"  -> universe = all genes in the sheet's gene column. Asks:
##                  "among the genes I analysed, is this column's set enriched
##                  for bucket X?" Correct when the gene column IS your analysis
##                  universe (e.g. all genes in the network). DEFAULT.
##      "syngo"  -> universe = all SynGO genes (within-synaptic-proteome test).
##      "genome" -> universe = all protein-coding genes (synaptic categories
##                  look enriched almost trivially if the sets are synaptic).
##      "custom" -> universe = CUSTOM_BACKGROUND (the set of genes that COULD
##                  have been hit / were tested). Best for GWAS/colocalisation.
##    => If the sheet only lists "interesting" genes rather than every analysed
##       gene, "sheet" tests a different (contrastive) question - switch to
##       "syngo" or "custom".
##  * PER-SOURCE BACKGROUND: SynGO buckets are tested against the SynGO-annotated
##    gene set (within-synaptic-proteome test); the GO glutamatergic category is
##    tested against the configured BACKGROUND_MODE universe. Universes (and N)
##    therefore differ by category source, so fold_enrichment/odds_ratio are not
##    directly comparable across sources and the BH families span both spaces.
##    The output carries 'source' (GO/SynGO) + 'source_detail' per category.
##  * GLUT_ONLY = TRUE restricts BOTH bucket definitions and the universe to
##    glutamatergic-synapse genes.
##  * Enrichment uses every bucket membership (not primary_bucket). Two FDRs are
##    reported: fdr_BH_within (across categories, per query set) and
##    fdr_BH_global (across all query_set x category tests).
##  * signalling_mediator is a residual/keyword bucket - low confidence. Read
##    its overlap genes before reporting.
## ============================================================================

## ---- 0. Packages -----------------------------------------------------------
cran_pkgs <- c("tidyverse", "readxl", "googlesheets4")
bioc_pkgs <- c("org.Hs.eg.db", "GO.db", "AnnotationDbi")

to_install <- cran_pkgs[!cran_pkgs %in% rownames(installed.packages())]
if (length(to_install)) install.packages(to_install)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
bioc_missing <- bioc_pkgs[!bioc_pkgs %in% rownames(installed.packages())]
if (length(bioc_missing)) BiocManager::install(bioc_missing, update = FALSE, ask = FALSE)

suppressPackageStartupMessages({
  library(tidyverse)        # dplyr, tidyr, stringr, purrr, tibble, readr, ...
  library(readxl)           # not attached by library(tidyverse)
  library(googlesheets4)    # not attached by library(tidyverse)
  library(org.Hs.eg.db); library(GO.db); library(AnnotationDbi)
})
select <- dplyr::select   # guard against AnnotationDbi masking dplyr::select

## ---- 1. CONFIG (edit these) ------------------------------------------------
SYNGO_ANNOTATION_XLSX <- "syngo_annotations.xlsx"   # unzipped SynGO release file
OUTPUT_TABLE_CSV      <- "glut_synapse_gene_categories.csv"
OUTPUT_ENRICH_CSV     <- "glut_synapse_enrichment.csv"
OUTPUT_OVERLAP_CSV    <- "glut_synapse_enrichment_overlap_genes.csv"
GLUT_ROOT_GO          <- "GO:0098978"               # glutamatergic synapse (CC)

## Results are also written back into the workbook as tabs (see WRITE_TO_SHEET).
WRITE_TO_SHEET        <- TRUE   # FALSE = local CSVs only, no write_sheet calls
OUTPUT_TABLE_SHEET    <- "glut synapse gene categories"
OUTPUT_ENRICH_SHEET   <- "glut synapse enrichment"
OUTPUT_OVERLAP_SHEET  <- "glut synapse enrichment overlap genes"

## --- Query gene sets: read from a Google Sheet ---
## You need EDIT access to this workbook (results are written back as tabs), so
## googlesheets4 will authenticate on first use (gs4_auth). Both reading and
## writing use that same authenticated session.
## Layout: one gene column + any number of TRUE/FALSE columns; each TRUE/FALSE
## column becomes a query set (genes flagged TRUE).
GOOGLE_SHEET_ID <- "<Google Sheet ID removed>"
SHEET_TAB       <- 1       # tab to read: name (e.g. "Sheet1") or 1-based position.
# NB: googlesheets4 targets by name/position, NOT gid.
GENE_COL        <- NULL    # gene-column name; NULL = auto-detect (gene/symbol)
QUERY_COLS      <- NULL    # restrict to these column names; NULL = all logical cols

BACKGROUND_MODE    <- "custom"       # "sheet" | "syngo" | "genome" | "custom"
CUSTOM_BACKGROUND  <- "ortholog_background_clean.txt"  # vector or file path; used iff mode=="custom"
GLUT_ONLY          <- FALSE          # restrict buckets + universe to glutamatergic
ENRICH_ALTERNATIVE <- "greater"      # "greater" (over-rep) | "two.sided"

## Bucket assignment by SynGO term-name regex (editable, inspectable).
bucket_name_patterns <- list(
  postsynaptic_receptor = "glutamate receptor|ionotropic|metabotropic|neurotransmitter receptor|AMPA|NMDA|kainate",
  presynaptic_release   = "active zone|synaptic vesicle|neurotransmitter (secretion|loading)|vesicle (docking|priming|fusion|exocytosis)|vesicle cycle|presynaptic .*release",
  scaffold              = "postsynaptic density|postsynaptic specialization|postsynaptic specialisation|structural constituent|scaffold|presynaptic cytomatrix",
  signalling_mediator   = "signal|signalling|signaling|second messenger|calcium ion|kinase|phosphat|G protein-coupled|GTPase|adenylate|modulation of chemical synaptic transmission"
)
curated_symbols <- list(
  postsynaptic_receptor = c(paste0("GRIN", c("1","2A","2B","2C","2D","3A","3B")),
                            paste0("GRIA", 1:4), paste0("GRIK", 1:5), paste0("GRM", 1:8)),
  scaffold              = c(paste0("DLG", 1:4), paste0("SHANK", 1:3), paste0("DLGAP", 1:5),
                            paste0("HOMER", 1:3), "GRIP1","GRIP2","PICK1","CASK","SYNGAP1"),
  presynaptic_release   = c("STX1A","STX1B","SNAP25","VAMP2","SYT1","SYT2",
                            paste0("RIMS", 1:4), "UNC13A","UNC13B","RAB3A","CPLX1","CPLX2",
                            "NSF","STXBP1", paste0("SLC17A", c("6","7","8")))
)
bucket_precedence <- c("postsynaptic_receptor","presynaptic_release","scaffold","signalling_mediator")

## --- Combined (conjunction) categories --------------------------------------
## Test enrichment for genes satisfying BOTH the GO glutamatergic flag AND a
## SynGO bucket (e.g. glutamatergic-synapse postsynaptic receptors). This is
## distinct from GLUT_ONLY (which restricts EVERY bucket + the universe to
## glutamatergic genes): here only the named combined category is a conjunction,
## while the standalone buckets are still tested in their normal spaces.
## Each entry:
##   bucket   = the SynGO bucket to intersect with the glutamatergic set
##   universe = denominator (drives the question being asked):
##     "syngo" -> SynGO proteome. Sits in the same space as the standalone bucket,
##                so the combined row is directly comparable and shows the
##                glutamatergic-specific increment. Most defensible. DEFAULT.
##     "glut"  -> glutamatergic-synapse genes only (conditional: "among
##                glutamatergic genes, is the query enriched for the bucket?").
##     "main"  -> the BACKGROUND_MODE universe (broadest; as per GO glutamatergic).
## Set to list() to disable. Combined rows carry source = "GO+SynGO".
COMBINED_CATEGORIES <- list(
  glutamatergic_x_postsynaptic_receptor = list(bucket = "postsynaptic_receptor",
                                               universe = "syngo")
)

## ---- helpers ---------------------------------------------------------------
load_gene_list <- function(x) {
  if (is.null(x)) return(character(0))
  if (length(x) == 1 && file.exists(x)) v <- readLines(x, warn = FALSE) else v <- as.character(x)
  v <- str_trim(v); v <- v[v != "" & !str_starts(v, "#")]; toupper(unique(v))
}
read_query_sheet <- function(id, tab = 1) {
  ## All columns read as character so the TRUE/FALSE parsing below stays robust.
  read_sheet(ss = id, sheet = tab, col_types = "c")
}
TRUE_TOKENS  <- c("TRUE","T","1","YES","Y","X")
FALSE_TOKENS <- c("FALSE","F","0","NO","N","")
as_logical_col <- function(x) {
  v <- toupper(str_trim(as.character(x)))
  out <- rep(NA, length(v)); out[v %in% TRUE_TOKENS] <- TRUE; out[v %in% FALSE_TOKENS] <- FALSE; out
}
is_logical_like <- function(x) {
  v <- toupper(str_trim(as.character(x))); v <- v[!is.na(v) & v != ""]
  length(v) > 0 && all(v %in% c(TRUE_TOKENS, FALSE_TOKENS))
}

## ---- 2. Glutamatergic gene universe from GO --------------------------------
glut_terms <- unique(c(GLUT_ROOT_GO, AnnotationDbi::get(GLUT_ROOT_GO, GOCCOFFSPRING)))
glut_terms <- glut_terms[!is.na(glut_terms)]
go2eg <- as.list(org.Hs.egGO2ALLEGS)
glut_entrez <- unique(unlist(go2eg[intersect(glut_terms, names(go2eg))]))
glut_entrez <- glut_entrez[!is.na(glut_entrez)]
glut_symbols <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = glut_entrez,
                                      column = "SYMBOL", keytype = "ENTREZID")
glut_symbols <- toupper(unique(na.omit(unname(glut_symbols))))
message(sprintf("Glutamatergic-synapse genes from GO: %d", length(glut_symbols)))

## ---- 3. Read & normalise SynGO annotations ---------------------------------
stopifnot(file.exists(SYNGO_ANNOTATION_XLSX))
syngo_raw <- read_excel(SYNGO_ANNOTATION_XLSX)
nm <- tolower(names(syngo_raw))
pick <- function(...) {
  for (p in c(...)) { hit <- which(str_detect(nm, p)); if (length(hit)) return(names(syngo_raw)[hit[1]]) }
  stop("No column matching: ", paste(c(...), collapse = " / "),
       "\nColumns: ", paste(names(syngo_raw), collapse = ", "))
}
col_symbol <- pick("hgnc_symbol", "^symbol$", "gene_symbol", "\\bsymbol\\b")
col_tname  <- pick("go.?term.?name", "term_name", "go_name", "\\bgo\\b.*name")
message("Resolved SynGO cols -> symbol: ", col_symbol, " | term_name: ", col_tname)

syngo <- syngo_raw %>%
  transmute(symbol    = toupper(str_trim(.data[[col_symbol]])),
            term_name = str_to_lower(str_trim(.data[[col_tname]]))) %>%
  filter(!is.na(symbol), symbol != "") %>% distinct()

## ---- 4. Location + buckets -> gene table -----------------------------------
location <- syngo %>%
  group_by(symbol) %>%
  summarise(presynaptic  = any(str_detect(term_name, "presynap"), na.rm = TRUE),
            postsynaptic = any(str_detect(term_name, "postsynap"), na.rm = TRUE),
            .groups = "drop")

bucket_from_terms <- map_dfr(names(bucket_name_patterns), function(b)
  syngo %>% filter(str_detect(term_name, regex(bucket_name_patterns[[b]], ignore_case = TRUE))) %>%
    distinct(symbol) %>% mutate(bucket = b))
bucket_from_symbols <- map_dfr(names(curated_symbols), function(b)
  tibble(symbol = toupper(curated_symbols[[b]]), bucket = b))
bucket_long <- bind_rows(bucket_from_terms, bucket_from_symbols) %>% distinct()

specific <- bucket_long %>%
  filter(bucket %in% c("postsynaptic_receptor","presynaptic_release","scaffold")) %>%
  pull(symbol) %>% unique()
bucket_long <- bucket_long %>% filter(!(bucket == "signalling_mediator" & symbol %in% specific))

bucket_wide <- bucket_long %>% mutate(val = TRUE) %>%
  pivot_wider(names_from = bucket, values_from = val, values_fill = FALSE)
for (b in bucket_precedence) if (!b %in% names(bucket_wide)) bucket_wide[[b]] <- FALSE
bucket_wide <- bucket_wide %>% rowwise() %>%
  mutate(n_buckets = sum(c_across(all_of(bucket_precedence))),
         primary_bucket = { hit <- bucket_precedence[as.logical(c_across(all_of(bucket_precedence)))]
         if (length(hit)) hit[1] else NA_character_ }) %>%
  ungroup() %>% mutate(multi_bucket = n_buckets > 1)

final <- tibble(symbol = sort(unique(c(syngo$symbol, bucket_long$symbol)))) %>%
  left_join(location, by = "symbol") %>%
  left_join(bucket_wide, by = "symbol") %>%
  mutate(glutamatergic = symbol %in% glut_symbols,
         across(c(presynaptic, postsynaptic, all_of(bucket_precedence), multi_bucket),
                ~ replace_na(.x, FALSE)),
         n_buckets = replace_na(n_buckets, 0L)) %>%
  relocate(symbol, glutamatergic, presynaptic, postsynaptic,
           primary_bucket, all_of(bucket_precedence), n_buckets, multi_bucket) %>%
  arrange(desc(glutamatergic), primary_bucket, symbol)
write_csv(final, OUTPUT_TABLE_CSV)
if (WRITE_TO_SHEET) write_sheet(final, ss = GOOGLE_SHEET_ID, sheet = OUTPUT_TABLE_SHEET)

## ---- 5. Read query sets from the Google Sheet ------------------------------
sheet <- read_query_sheet(GOOGLE_SHEET_ID, SHEET_TAB)
message(sprintf("Sheet read: %d rows x %d cols (%s)",
                nrow(sheet), ncol(sheet), paste(names(sheet), collapse = ", ")))

gene_col <- GENE_COL
if (is.null(gene_col)) {
  cand <- names(sheet)[str_detect(tolower(names(sheet)), "gene|symbol")]
  gene_col <- if (length(cand)) cand[1] else names(sheet)[1]
}
message("Gene column: ", gene_col)
genes_vec  <- toupper(str_trim(as.character(sheet[[gene_col]])))
sheet_genes <- unique(genes_vec[!is.na(genes_vec) & genes_vec != ""])

cand_cols <- setdiff(names(sheet), gene_col)
if (!is.null(QUERY_COLS)) cand_cols <- intersect(QUERY_COLS, cand_cols)
logical_cols <- cand_cols[map_lgl(sheet[cand_cols], is_logical_like)]
skipped <- setdiff(cand_cols, logical_cols)
if (length(skipped)) message("Skipping non-TRUE/FALSE columns: ", paste(skipped, collapse = ", "))
if (!length(logical_cols)) stop("No TRUE/FALSE query columns found in the sheet.")

query_sets <- logical_cols %>%
  set_names() %>%
  map(function(cl) {
    flag <- as_logical_col(sheet[[cl]])
    g <- unique(genes_vec[which(flag %in% TRUE)])
    g[!is.na(g) & g != ""]
  })
message("Query sets (", length(query_sets), "): ",
        paste(sprintf("%s[%d]", names(query_sets), lengths(query_sets)), collapse = ", "))

## ---- 6. Define per-source background universes -----------------------------
## SynGO buckets -> SynGO-annotated genes; GO glutamatergic -> BACKGROUND_MODE
## universe. Categories thus carry their own universe + N (see methods note).
main_universe <- switch(BACKGROUND_MODE,
                        sheet  = sheet_genes,
                        syngo  = unique(syngo$symbol),
                        genome = toupper(unique(na.omit(AnnotationDbi::keys(org.Hs.eg.db, keytype = "SYMBOL")))),
                        custom = { u <- load_gene_list(CUSTOM_BACKGROUND)
                        if (!length(u)) stop("BACKGROUND_MODE='custom' but CUSTOM_BACKGROUND is empty."); u },
                        stop("Unknown BACKGROUND_MODE: ", BACKGROUND_MODE)
)
syngo_universe <- unique(syngo$symbol)          # all genes annotated in SynGO
if (GLUT_ONLY) {
  main_universe  <- intersect(main_universe,  glut_symbols)
  syngo_universe <- intersect(syngo_universe, glut_symbols)
}
main_universe  <- unique(main_universe)
syngo_universe <- unique(syngo_universe)
message(sprintf("Main background = '%s' -> %d genes%s",
                BACKGROUND_MODE, length(main_universe),
                if (GLUT_ONLY) " (glutamatergic-only)" else ""))
message(sprintf("SynGO background -> %d genes (all SynGO-annotated)%s",
                length(syngo_universe), if (GLUT_ONLY) " (glutamatergic-only)" else ""))

## Per-category spec: tested gene set, its own universe, source label + detail.
## SynGO bucket genes are intersected with the SynGO universe, so curated symbols
## absent from SynGO are dropped from both category and background (consistent).
cat_specs <- list()
for (b in bucket_precedence) {
  cat_specs[[b]] <- list(
    genes    = intersect(final$symbol[final[[b]]], syngo_universe),
    universe = syngo_universe,
    source   = "SynGO",
    detail   = "SynGO term-name regex + curated synaptic symbols"
  )
}
if (!GLUT_ONLY) {
  cat_specs[["glutamatergic"]] <- list(
    genes    = intersect(glut_symbols, main_universe),
    universe = main_universe,
    source   = "GO",
    detail   = paste0(GLUT_ROOT_GO, " (glutamatergic synapse, CC) + GOCCOFFSPRING")
  )
}

## --- Combined conjunction categories (GO glutamatergic AND a SynGO bucket) ---
## Numerator = glutamatergic genes AND bucket members (bucket defined in SynGO
## space, so a glut gene absent from SynGO counts as a non-member, not the
## category). Denominator per COMBINED_CATEGORIES$universe. Under GLUT_ONLY the
## chosen universe is already glut-restricted upstream, so it collapses sensibly.
combined_universe_of <- function(key) {
  u <- switch(key,
              syngo = syngo_universe,
              glut  = glut_symbols,
              main  = main_universe,
              stop("Unknown combined-category universe: '", key, "'"))
  if (GLUT_ONLY) u <- intersect(u, glut_symbols)
  unique(u)
}
for (nm in names(COMBINED_CATEGORIES)) {
  spec_in <- COMBINED_CATEGORIES[[nm]]
  b <- spec_in$bucket
  if (!b %in% bucket_precedence) stop("COMBINED_CATEGORIES['", nm, "']: unknown bucket '", b, "'")
  uni          <- combined_universe_of(spec_in$universe)
  bucket_genes <- intersect(final$symbol[final[[b]]], syngo_universe)  # bucket lives in SynGO space
  genes        <- Reduce(intersect, list(glut_symbols, bucket_genes, uni))
  cat_specs[[nm]] <- list(
    genes    = genes,
    universe = uni,
    source   = "GO+SynGO",
    detail   = sprintf("glutamatergic (%s + GOCCOFFSPRING) AND SynGO bucket '%s'; universe = %s",
                       GLUT_ROOT_GO, b, spec_in$universe)
  )
  message(sprintf("Combined category '%s': %d genes in a universe of %d (universe = %s)",
                  nm, length(genes), length(uni), spec_in$universe))
}

## ---- 7. Over-representation analysis (per query set) -----------------------
## Each category is tested against ITS OWN universe (spec$universe), so N and the
## in-universe query size (query_k) vary by category/source.
run_ora <- function(query_all, label) {
  imap_dfr(cat_specs, function(spec, cat) {
    uni     <- spec$universe
    N       <- length(uni)
    query   <- intersect(query_all, uni)
    k       <- length(query)
    genes   <- spec$genes
    m       <- length(genes)
    overlap <- intersect(query, genes); a <- length(overlap)
    mat <- matrix(c(a, k - a, m - a, N - k - m + a), nrow = 2, byrow = TRUE)
    ft  <- fisher.test(mat, alternative = ENRICH_ALTERNATIVE)
    tibble(query_set = label, category = cat,
           source = spec$source, source_detail = spec$detail,
           universe_N = N, query_supplied = length(query_all), query_k = k,
           query_dropped = length(setdiff(query_all, uni)),
           bucket_m = m, overlap = a,
           expected = round(k * m / N, 2),
           fold_enrichment = ifelse(m == 0 | k == 0, NA, round((a / k) / (m / N), 2)),
           odds_ratio = round(unname(ft$estimate), 2),
           p_value = ft$p.value, overlap_genes = paste(sort(overlap), collapse = ";"))
  })
}
message("\nRunning enrichment per query set (per-source backgrounds):")
enrich <- imap_dfr(query_sets, run_ora) %>%
  group_by(query_set) %>% mutate(fdr_BH_within = p.adjust(p_value, "BH")) %>% ungroup() %>%
  mutate(fdr_BH_global = p.adjust(p_value, "BH")) %>%
  arrange(query_set, p_value) %>%
  relocate(query_set, category, source, source_detail, overlap, expected,
           fold_enrichment, odds_ratio, p_value, fdr_BH_within, fdr_BH_global)

write_csv(select(enrich, -overlap_genes), OUTPUT_ENRICH_CSV)
write_csv(select(enrich, query_set, category, source, overlap, overlap_genes), OUTPUT_OVERLAP_CSV)
if (WRITE_TO_SHEET) {
  write_sheet(select(enrich, -overlap_genes),
              ss = GOOGLE_SHEET_ID, sheet = OUTPUT_ENRICH_SHEET)
  write_sheet(select(enrich, query_set, category, source, overlap, overlap_genes),
              ss = GOOGLE_SHEET_ID, sheet = OUTPUT_OVERLAP_SHEET)
}

## ---- 8. Report -------------------------------------------------------------
message("\n==== ENRICHMENT (SynGO buckets vs SynGO background; ",
        "GO glutamatergic vs '", BACKGROUND_MODE, "'; ",
        "combined GO+SynGO categories per their configured universe",
        if (GLUT_ONLY) ", glutamatergic-only" else "", ") ====")
print(as.data.frame(select(enrich, query_set, category, source, overlap, expected,
                           fold_enrichment, odds_ratio, p_value,
                           fdr_BH_within, fdr_BH_global)),
      row.names = FALSE, digits = 3)
message("\nWrote CSVs:\n  ", OUTPUT_TABLE_CSV, "\n  ", OUTPUT_ENRICH_CSV, "\n  ", OUTPUT_OVERLAP_CSV)
if (WRITE_TO_SHEET)
  message("Wrote sheet tabs:\n  ", OUTPUT_TABLE_SHEET, "\n  ", OUTPUT_ENRICH_SHEET,
          "\n  ", OUTPUT_OVERLAP_SHEET)