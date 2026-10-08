# ============================================================
# Query GWAS Catalog hits for genes and annotate traits with ChatGPT
# v3: temperament annotations added; sequential lookup with 10s timeout
# ============================================================

library(httr)
library(jsonlite)
library(dplyr)
library(purrr)
library(tibble)
library(stringr)
library(readr)

#========================
# Settings
#========================

GENES <- c(
  "ACSF2","ADCY2","ADGRA1","ADGRB3","ADGRD1","ADGRF3","ADGRG1","ADGRG2","ADGRG7",
  "ADRA2C","AEBP2","AJAP1","AKAP6","ANKIB1","ANTXR2","ARGLU1","ARHGAP15","ARHGAP33",
  "ARIH1","ARIH2","ATL2","ATP1A2","ATP5MC1","AURKB","B3GAT2","BIK","BNC1","BSN",
  "BTN3A2","BTN3A3","C20orf144","C2CD3","C3orf62","CA9","CACNA1H","CACNB2","CACNG3",
  "CADM2","CAPN7","CCDC71","CCHCR1","CD34","CDH10","CDH12","CDH18","CDH2","CDH8",
  "CDH9","CELF4","CELF5","CELSR3","CHL1","CHST9","CLP1","CNGB1","CNTN1","CNTN6",
  "CNTNAP1","CNTNAP4","COL6A5","CRTC1","CTNND1","CYP1B1","DAAM2","DACH1","DALRD3",
  "DAOA","DCLK1","DGCR6L","DGKI","DHX35","DLG2","DLGAP1","DPP10","DUOX2","EDEM1",
  "ELAVL2","ELAVL3","ELK3","EMC2","ENC1","EPHA4","EVI2B","F2RL1","FAM155A","FHOD1",
  "FLOT1","FMNL3","GABRA1","GABRA2","GABRA5","GABRA6","GABRB2","GABRB3","GABRG1",
  "GALK2","GAP43","GAS7","GATAD2A","GK2","GMDS","GNPDA2","GPR137","GRIA1","GRIA2",
  "GRIA3","GRIA4","GRIK1","GRIK2","GRIK3","GRIN2A","GRIN2B","GRM3","GRM4","GRM5",
  "GRM7","H4C11","HCN1","HCN2","HLA-B","HLA-DMA","HS3ST2","HSPB6","ID2","IER3",
  "IGFLR1","IGSF11","IL12RB2","IQCC","IQSEC3","IVNS1ABP","KANSL2","KATNB1","KATNIP",
  "KBTBD12","KCNB1","KCNC1","KCND2","KCND3","KCNE4","KCNIP1","KCNV1","KLHDC8B",
  "KLHL1","KLHL22","KLHL24","KLHL26","KLHL31","KLHL33","KLHL35","KLHL38","KLHL40",
  "KLHL41","KLHL5","KLHL6","KLHL7","KMT2B","L1CAM","LAMB2","LIN37","LINGO1","LRFN2",
  "LRFN3","LRFN5","LRPAP1","LRRC4C","LRRN1","LRRTM1","LRRTM2","LSAMP","LY86","MAFB",
  "MAIP1","MEF2C","MLN","MTSS2","NCAM2","NCKIPSD","NECAB1","NECAB3","NES","NFASC",
  "NLGN1","NMUR2","NOVA1","NPHP4","NPHS1","NR2C2","NRCAM","NRGN","NRM","NTRK2","NTRK3",
  "NYNRIN","OGFRL1","OMG","OPCML","P4HTM","PABPC1L","PAK3","PAWR","PBX1","PCDH10",
  "PCDH19","PCSK1N","PDE1B","PDE3A","PDE4D","PGBD1","PLAAT4","PLXNA2","PLXNA4","PPA2",
  "PPFIA2","PRODH2","PROSER3","PRPH2","PSENEN","PWWP2A","QRICH1","RABEPK","RASGRF1",
  "RBFOX1","RBL1","RHOB","RNASEH2A","RNF144A","RORB","RREB1","RTN4R","SAMD9L","SCG2",
  "SCN2A","SCN2B","SEMA6B","SEMA6C","SEMA6D","SEPTIN3","SERBP1","SERPING1","SH3GL3",
  "SIGIRR","SLC25A17","SLC4A10","SLC4A3","SLC6A11","SLC6A15","SLIT2","SPDYC","SPMIP2",
  "SPSB1","SPSB4","SRGAP3","STXBP1","SYT4","TENT4A","TEX9","TMEM107","TMPRSS15","TMX2",
  "TRHR","TRIM27","TUBB","U2AF1L4","UBA7","UBE2L6","UBE2Z","UBR2","UNC80","USP4","UTP20",
  "VARS2","VSNL1","WDR6","XPNPEP3","YPEL4","YWHAB","ZDHHC5","ZFHX3","ZKSCAN4","ZNF396",
  "ZNF804B","ZNF839","ZSCAN16","ZSCAN23","ZWINT","CBLL2","EP300","ZNF804A","DGAT2L6",
  "CLEC2L","RGCC","UBR3","ZNF362","GPR179","ERAP1","SECTM1","TREX2","CAMLG","SLC9A5",
  "ZCCHC14","DGCR6","SPINDOC","FAXC","ZNF408","ARPIN","SKIDA1","FAM228B")

GWS_THRESHOLD <- 5e-8
BASE_URL      <- "https://www.ebi.ac.uk/gwas/rest/api"
CACHE_DIR     <- "gwas_cache"
dir.create(CACHE_DIR, showWarnings = FALSE)

gene_assoc_cache_file   <- file.path(CACHE_DIR, "all_gene_assocs.csv")
gene_status_file        <- file.path(CACHE_DIR, "gene_query_status.csv")
gws_cache_file          <- file.path(CACHE_DIR, "gws_hits.csv")
study_cache_file        <- file.path(CACHE_DIR, "study_info.csv")
assoc_study_cache_file  <- file.path(CACHE_DIR, "assoc_to_study.csv")
# Reuse v2 trait annotation cache — same schema
trait_annotation_cache_file <- file.path(CACHE_DIR, "trait_annotations_chatgpt.v2.csv")

# Output file read by parse_gwas_catalog_hits.R
out_file <- "gwas_catalog_hits.csv"

chatgpt_model    <- "gpt-4o-mini"
trait_batch_size <- 25
checkpoint_every <- 100   # association-study lookups per cache write

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b

#========================
# General helpers
#========================

empty_assoc_tbl <- function() {
  tibble(
    variant_id = character(), association_id = character(), study_id = character(),
    study_url = character(), pvalue = numeric(), pvalue_mantissa = numeric(),
    pvalue_exponent = numeric(), or_per_copy = numeric(), beta_num = numeric(),
    beta_unit = character(), queried_gene = character()
  )
}

standardize_assoc_cols <- function(x) {
  if (is.null(x) || nrow(x) == 0) return(empty_assoc_tbl())
  for (nm in names(empty_assoc_tbl())) if (!nm %in% names(x)) x[[nm]] <- NA
  x %>%
    mutate(
      variant_id     = as.character(variant_id),
      association_id = as.character(association_id),
      study_id       = as.character(study_id),
      study_url      = as.character(study_url),
      pvalue         = as.numeric(pvalue),
      pvalue_mantissa = as.numeric(pvalue_mantissa),
      pvalue_exponent = as.numeric(pvalue_exponent),
      or_per_copy    = as.numeric(or_per_copy),
      beta_num       = as.numeric(beta_num),
      beta_unit      = as.character(beta_unit),
      queried_gene   = as.character(queried_gene)
    ) %>%
    select(all_of(names(empty_assoc_tbl())))
}

write_gene_status <- function(all_assocs, genes = GENES) {
  per_gene_cached_genes <- str_remove(
    basename(list.files(CACHE_DIR, pattern = "_assocs\\.csv$", full.names = TRUE)),
    "_assocs\\.csv$"
  )
  gene_query_status <- tibble(queried_gene = genes) %>%
    mutate(
      per_gene_cache_exists = queried_gene %in% per_gene_cached_genes,
      in_all_assocs         = queried_gene %in% unique(all_assocs$queried_gene),
      n_all_associations    = map_int(queried_gene,
                                      function(g) sum(all_assocs$queried_gene == g, na.rm = TRUE))
    )
  write_csv(gene_query_status, gene_status_file)
  gene_query_status
}

#========================
# GWAS Catalog helpers
#========================

get_snps_for_gene <- function(gene) {
  url  <- paste0(BASE_URL, "/singleNucleotidePolymorphisms/search/findByGene",
                 "?geneName=", gene, "&size=1000")
  resp <- GET(url, timeout(60))
  if (status_code(resp) != 200) return(NULL)
  json <- content(resp, as = "text", encoding = "UTF-8") %>% fromJSON(flatten = TRUE)
  snps <- json$`_embedded`$singleNucleotidePolymorphisms
  if (is.null(snps) || !is.data.frame(snps) || nrow(snps) == 0) return(NULL)
  snps %>%
    select(rsId, assoc_url = `_links.associationsBySnpSummary.href`) %>%
    distinct(rsId, .keep_all = TRUE)
}

get_assocs_for_snp <- function(rsid, assoc_url) {
  resp  <- GET(assoc_url, timeout(60))
  if (status_code(resp) != 200) return(NULL)
  json  <- content(resp, as = "text", encoding = "UTF-8") %>% fromJSON(flatten = TRUE)
  assocs <- json$`_embedded`$associations
  if (is.null(assocs) || !is.data.frame(assocs) || nrow(assocs) == 0) return(NULL)
  tryCatch({
    tibble(
      variant_id      = rsid,
      association_id  = tryCatch(basename(assocs$`_links.self.href`),          error = function(e) rep(NA_character_, nrow(assocs))),
      study_id        = NA_character_,
      study_url       = tryCatch(as.character(assocs$`_links.study.href`),      error = function(e) rep(NA_character_, nrow(assocs))),
      pvalue          = tryCatch(as.numeric(assocs$pvalue),                     error = function(e) rep(NA_real_,      nrow(assocs))),
      pvalue_mantissa = tryCatch(as.numeric(assocs$pvalueMantissa),             error = function(e) rep(NA_real_,      nrow(assocs))),
      pvalue_exponent = tryCatch(as.numeric(assocs$pvalueExponent),             error = function(e) rep(NA_real_,      nrow(assocs))),
      or_per_copy     = tryCatch(as.numeric(assocs$orPerCopyNum),               error = function(e) rep(NA_real_,      nrow(assocs))),
      beta_num        = tryCatch(as.numeric(assocs$betaNum),                    error = function(e) rep(NA_real_,      nrow(assocs))),
      beta_unit       = tryCatch(as.character(assocs$betaUnit),                 error = function(e) rep(NA_character_, nrow(assocs)))
    )
  }, error = function(e) {
    warning("Error parsing associations for ", rsid, ": ", conditionMessage(e))
    NULL
  })
}

get_study_info <- function(study_id) {
  url  <- paste0(BASE_URL, "/studies/", study_id)
  resp <- GET(url, timeout(60))
  if (status_code(resp) != 200) return(NULL)
  json <- content(resp, as = "text", encoding = "UTF-8") %>% fromJSON(flatten = FALSE)
  tibble(
    study_id       = study_id,
    reported_trait = tryCatch(json$diseaseTrait$trait,              error = function(e) NA_character_),
    pubmed_id      = tryCatch(json$publicationInfo$pubmedId,        error = function(e) NA_character_),
    title          = tryCatch(json$publicationInfo$title,           error = function(e) NA_character_),
    author         = tryCatch(json$publicationInfo$author$fullname, error = function(e) NA_character_),
    pub_date       = tryCatch(json$publicationInfo$publicationDate, error = function(e) NA_character_),
    sample_size    = tryCatch(json$initialSampleSize,               error = function(e) NA_character_)
  ) %>% standardize_study_cols()
}

study_info_cols <- c("study_id","reported_trait","pubmed_id","title","author","pub_date","sample_size")

empty_study_tbl <- function() {
  tibble(study_id=character(), reported_trait=character(), pubmed_id=character(),
         title=character(), author=character(), pub_date=character(), sample_size=character())
}

standardize_study_cols <- function(x) {
  if (is.null(x) || nrow(x) == 0) return(empty_study_tbl())
  for (nm in study_info_cols) if (!nm %in% names(x)) x[[nm]] <- NA_character_
  x %>% mutate(across(all_of(study_info_cols), as.character)) %>% select(all_of(study_info_cols))
}

# Note: study_url is /associations/{id}/study — no GCST in URL, API call required.
get_study_id_for_assoc <- function(association_id) {
  r <- GET(paste0(BASE_URL, "/associations/", association_id, "/study"), timeout(10))
  if (status_code(r) != 200) return(NA_character_)
  j <- content(r, as = "text", encoding = "UTF-8") %>% fromJSON(flatten = TRUE)
  j$accessionId %||% str_extract(j$`_links.self.href` %||% "", "GCST[0-9]+")
}

#========================
# ChatGPT helpers
#========================

strip_json_fences <- function(x) {
  x %>% str_trim() %>%
    str_remove("^```json\\s*") %>% str_remove("^```\\s*") %>% str_remove("\\s*```$")
}

extract_response_text <- function(result) {
  if (!is.null(result$output_text)) {
    out <- paste(result$output_text, collapse = "\n")
    if (nchar(trimws(out)) > 0) return(out)
  }
  if (!is.null(result$output) && is.data.frame(result$output) && "content" %in% names(result$output)) {
    for (i in seq_len(nrow(result$output))) {
      content_i <- result$output$content[[i]]
      if (is.null(content_i)) next
      if (is.data.frame(content_i) && "text" %in% names(content_i)) {
        out <- paste(content_i$text, collapse = "\n")
        if (nchar(trimws(out)) > 0) return(out)
      }
      if (is.list(content_i)) {
        text_blocks <- unlist(lapply(content_i, function(y) {
          if (is.list(y) && !is.null(y$text)) return(y$text)
          if (is.data.frame(y) && "text" %in% names(y)) return(y$text)
          NULL
        }), use.names = FALSE)
        out <- paste(text_blocks, collapse = "\n")
        if (nchar(trimws(out)) > 0) return(out)
      }
    }
  }
  ""
}

call_chatgpt <- function(system_prompt, user_prompt, model = chatgpt_model, max_tries = 3) {
  if (Sys.getenv("OPENAI_API_KEY") == "")
    stop("OPENAI_API_KEY is not set. Run: Sys.setenv(OPENAI_API_KEY = 'your_key_here')")
  for (try_i in seq_len(max_tries)) {
    message("    ChatGPT API attempt ", try_i, "/", max_tries)
    body_list <- list(
      model  = model,
      input  = list(list(role = "system", content = system_prompt),
                    list(role = "user",   content = user_prompt)),
      max_output_tokens = 8000
    )
    resp <- POST(
      url     = "https://api.openai.com/v1/responses",
      add_headers("Content-Type" = "application/json",
                  "Authorization" = paste("Bearer", Sys.getenv("OPENAI_API_KEY"))),
      body    = toJSON(body_list, auto_unbox = TRUE),
      encode  = "raw",
      timeout(120)
    )
    resp_text <- content(resp, as = "text", encoding = "UTF-8")
    if (status_code(resp) == 200) {
      result   <- fromJSON(resp_text, flatten = FALSE)
      out_text <- extract_response_text(result)
      if (!is.null(out_text) && nchar(trimws(out_text)) > 0) return(out_text)
      message("    Empty parsed model response. Full API response:\n", resp_text)
    } else {
      message("    OpenAI API error ", status_code(resp), ": ", resp_text)
    }
    if (try_i < max_tries) Sys.sleep(5 * try_i)
  }
  stop("OpenAI API call failed or returned empty text after ", max_tries, " attempts.")
}

coerce_logical_clean <- function(x) {
  x_chr <- tolower(as.character(x))
  case_when(x_chr == "true" ~ TRUE, x_chr == "false" ~ FALSE, TRUE ~ NA)
}

clean_null <- function(x) { x <- as.character(x); na_if(na_if(na_if(x, "NULL"), "null"), "") }

annotate_trait_batch_with_chatgpt <- function(traits, batch_i, n_batches) {
  traits <- traits[!is.na(traits) & traits != ""]
  if (length(traits) == 0) {
    return(tibble(trait=character(), is_psychiatric=logical(), disorder=character(),
                  is_temperament=logical(), temperament_category=character()))
  }
  trait_list_str <- paste(seq_along(traits), traits, sep = ". ", collapse = "\n")

  system_prompt <- paste0(
    "You are annotating GWAS reported traits.\n",
    "Return ONLY a JSON array. No markdown. No prose.\n",
    "Each element must have exactly these keys: ",
    "trait, is_psychiatric, disorder, is_temperament, temperament_category."
  )
  user_prompt <- paste0(
    "For each trait below, return a JSON array where each element has:\n",
    "  - \"trait\": the exact trait string\n",
    "  - \"is_psychiatric\": true if the trait is a psychiatric, neurological, ",
    "or behavioural disorder or symptom; otherwise false\n",
    "  - \"disorder\": one of \"depression\", \"schizophrenia\", \"OCD\", ",
    "\"other_psychiatric\", or null\n",
    "  - \"is_temperament\": true if the trait measures personality, temperament, or ",
    "a related stable psychological disposition (see categories below); otherwise false\n",
    "  - \"temperament_category\": one of the following or null:\n",
    "      \"Big Five\" — openness, conscientiousness, extraversion/introversion, ",
    "agreeableness, neuroticism, or facets thereof\n",
    "      \"temperament\" — harm avoidance, novelty seeking, reward dependence, ",
    "persistence, behavioural inhibition/activation, impulsivity, irritability, ",
    "sensation seeking, risk tolerance, emotional reactivity, or similar stable traits\n",
    "      \"personality_other\" — other personality measures not captured above ",
    "(e.g. self-directedness, cooperativeness, well-being, life satisfaction)\n",
    "      null — if is_temperament is false\n\n",
    "Rules:\n",
    "- disorder = \"depression\" if it relates to major depression, depressive episodes, or MDD\n",
    "- disorder = \"schizophrenia\" if it relates to schizophrenia or schizoaffective disorder\n",
    "- disorder = \"OCD\" if it relates to obsessive-compulsive disorder\n",
    "- disorder = \"other_psychiatric\" if psychiatric but not one of the three above\n",
    "- disorder = null if not psychiatric\n",
    "- is_temperament and is_psychiatric can both be true\n",
    "- Return one JSON object for every input trait\n\n",
    "Traits:\n", trait_list_str
  )

  message("\nChatGPT trait batch ", batch_i, "/", n_batches, " | traits: ", length(traits))
  gpt_text <- call_chatgpt(system_prompt, user_prompt, model = chatgpt_model)
  json_str  <- strip_json_fences(gpt_text)

  annotations <- tryCatch(fromJSON(json_str, flatten = TRUE), error = function(e) {
    warning("Failed to parse ChatGPT JSON for batch ", batch_i, ": ",
            conditionMessage(e), "\nRaw response:\n", json_str)
    NULL
  })
  if (is.null(annotations)) stop("ChatGPT response could not be parsed for batch ", batch_i)

  annotations <- as_tibble(annotations)
  required_cols <- c("trait","is_psychiatric","disorder","is_temperament","temperament_category")
  missing_cols  <- setdiff(required_cols, names(annotations))
  if (length(missing_cols) > 0)
    stop("Missing expected columns in ChatGPT response: ", paste(missing_cols, collapse = ", "))

  annotations <- annotations %>%
    transmute(
      trait                = as.character(trait),
      is_psychiatric       = coerce_logical_clean(is_psychiatric),
      disorder             = clean_null(disorder),
      is_temperament       = coerce_logical_clean(is_temperament),
      temperament_category = clean_null(temperament_category)
    )

  missing_traits <- setdiff(traits, annotations$trait)
  if (length(missing_traits) > 0)
    warning("Some traits not returned by ChatGPT:\n", paste(missing_traits, collapse = "\n"))

  annotations
}

annotate_traits_with_chatgpt_cached <- function(traits) {
  unique_traits <- unique(traits[!is.na(traits) & traits != ""])
  empty_cache   <- tibble(trait=character(), is_psychiatric=logical(), disorder=character(),
                          is_temperament=logical(), temperament_category=character())
  if (length(unique_traits) == 0) return(empty_cache)

  if (file.exists(trait_annotation_cache_file)) {
    cached <- read_csv(trait_annotation_cache_file, show_col_types = FALSE) %>%
      mutate(
        trait                = as.character(trait),
        is_psychiatric       = coerce_logical_clean(is_psychiatric),
        disorder             = clean_null(disorder),
        is_temperament       = if ("is_temperament" %in% names(.)) coerce_logical_clean(is_temperament) else NA,
        temperament_category = if ("temperament_category" %in% names(.)) clean_null(temperament_category) else NA_character_
      )
  } else {
    cached <- empty_cache
  }

  cached <- cached %>% filter(!is.na(trait), trait != "") %>% distinct(trait, .keep_all = TRUE)
  traits_needed <- setdiff(unique_traits, cached$trait)
  message("Trait annotations cached: ", nrow(cached))
  message("Trait annotations needed: ", length(traits_needed))

  if (length(traits_needed) > 0) {
    batches <- split(traits_needed, ceiling(seq_along(traits_needed) / trait_batch_size))
    for (batch_i in seq_along(batches)) {
      new_annotations <- annotate_trait_batch_with_chatgpt(batches[[batch_i]], batch_i, length(batches))
      cached <- bind_rows(cached, new_annotations) %>% distinct(trait, .keep_all = TRUE)
      write_csv(cached, trait_annotation_cache_file)
      message("  Saved annotation checkpoint: ", nrow(cached), " traits annotated")
    }
  }
  cached %>% filter(trait %in% unique_traits)
}

#========================
# Step 1: Query GWAS Catalog associations, cached per gene
#========================

message("Querying GWAS Catalog for ", length(GENES), " genes...")

if (file.exists(gene_assoc_cache_file)) {
  message("Loading all gene associations cache: ", gene_assoc_cache_file)
  all_assocs <- read_csv(gene_assoc_cache_file, show_col_types = FALSE) %>% standardize_assoc_cols()
} else {
  all_assocs <- empty_assoc_tbl()
}

per_gene_cached_genes <- str_remove(
  basename(list.files(CACHE_DIR, pattern = "_assocs\\.csv$", full.names = TRUE)), "_assocs\\.csv$")
cached_genes  <- union(unique(all_assocs$queried_gene), per_gene_cached_genes)
genes_needed  <- setdiff(GENES, cached_genes)

message("Genes already cached: ", length(cached_genes))
message("Genes still needing query: ", length(genes_needed))

for (gene in genes_needed) {
  cache_file_gene <- file.path(CACHE_DIR, paste0(gene, "_assocs.csv"))
  message(" -> ", gene, appendLF = FALSE)
  snps <- get_snps_for_gene(gene)
  if (is.null(snps)) {
    message(" [no SNPs]"); write_csv(empty_assoc_tbl(), cache_file_gene); write_gene_status(all_assocs); next
  }
  gene_assocs <- map_dfr(seq_len(nrow(snps)), function(i) {
    df <- tryCatch(get_assocs_for_snp(snps$rsId[i], snps$assoc_url[i]), error = function(e) NULL)
    if (!is.null(df) && nrow(df) > 0) df$queried_gene <- gene
    df
  }) %>% standardize_assoc_cols()
  if (nrow(gene_assocs) == 0) {
    message(" [no associations]"); write_csv(empty_assoc_tbl(), cache_file_gene); write_gene_status(all_assocs); next
  }
  write_csv(gene_assocs, cache_file_gene)
  all_assocs <- bind_rows(standardize_assoc_cols(all_assocs), standardize_assoc_cols(gene_assocs)) %>%
    distinct(queried_gene, variant_id, association_id, study_url, pvalue, .keep_all = TRUE)
  write_csv(all_assocs, gene_assoc_cache_file)
  write_gene_status(all_assocs)
  message(" [", nrow(gene_assocs), " associations]")
}

for (gene in GENES) {
  cache_file_gene <- file.path(CACHE_DIR, paste0(gene, "_assocs.csv"))
  if (!file.exists(cache_file_gene) || gene %in% unique(all_assocs$queried_gene)) next
  gene_assocs <- read_csv(cache_file_gene, show_col_types = FALSE) %>% standardize_assoc_cols()
  if (nrow(gene_assocs) > 0) {
    gene_assocs$queried_gene <- gene
    all_assocs <- bind_rows(all_assocs, gene_assocs) %>% standardize_assoc_cols() %>%
      distinct(queried_gene, variant_id, association_id, study_url, pvalue, .keep_all = TRUE)
    write_csv(all_assocs, gene_assoc_cache_file)
  }
}

gene_query_status <- write_gene_status(all_assocs)
message("Final all_assocs genes represented: ", n_distinct(all_assocs$queried_gene))
message("Genes with per-gene cache files: ",    sum(gene_query_status$per_gene_cache_exists))
message("Genes with association rows: ",        sum(gene_query_status$n_all_associations > 0))

#========================
# Step 2: Filter to genome-wide significant hits
#========================

gws <- all_assocs %>%
  filter(!is.na(pvalue), pvalue <= GWS_THRESHOLD, !is.na(association_id), association_id != "") %>%
  mutate(study_id = NA_character_)

message("Genome-wide significant associations: ", nrow(gws))

#========================
# Step 2b: Recover study IDs via API lookup, cached per association
#
# study_url is /associations/{id}/study — no GCST in the URL itself,
# so an API call is needed. This is a one-time cost; results are cached.
#========================

if (file.exists(assoc_study_cache_file)) {
  assoc_to_study <- read_csv(assoc_study_cache_file, show_col_types = FALSE) %>%
    mutate(association_id = as.character(association_id), study_id = as.character(study_id))
} else {
  assoc_to_study <- tibble(association_id = character(), study_id = character())
}

assoc_needed <- setdiff(gws %>% distinct(association_id) %>% pull(), assoc_to_study$association_id)

message("Association-study cache entries: ",     nrow(assoc_to_study))
message("Association-study lookups still needed: ", length(assoc_needed))

if (length(assoc_needed) > 0) {
  message("Running ", length(assoc_needed), " association-study lookups...")
  for (i in seq_along(assoc_needed)) {
    aid <- assoc_needed[i]
    new_row <- tibble(
      association_id = as.character(aid),
      study_id       = tryCatch(get_study_id_for_assoc(aid), error = function(e) NA_character_)
    )
    assoc_to_study <- bind_rows(assoc_to_study, new_row) %>%
      mutate(association_id = as.character(association_id), study_id = as.character(study_id)) %>%
      distinct(association_id, .keep_all = TRUE)
    if (i %% checkpoint_every == 0 || i == length(assoc_needed)) {
      write_csv(assoc_to_study, assoc_study_cache_file)
      message("  ", i, "/", length(assoc_needed), " lookups done — cache saved")
    }
  }
  message("Association-study lookups complete.")
}

gws <- gws %>%
  left_join(assoc_to_study, by = "association_id", suffix = c("", ".fallback")) %>%
  mutate(study_id = coalesce(study_id, study_id.fallback)) %>%
  select(-any_of("study_id.fallback")) %>%
  filter(!is.na(study_id), study_id != "", grepl("^GCST", study_id))

write_csv(gws, gws_cache_file)
message("GWS associations with study IDs: ", nrow(gws))

#========================
# Step 3: Fetch study info, cached per study
#========================

unique_study_ids <- unique(gws$study_id)
message("Fetching info for ", length(unique_study_ids), " unique studies...")

if (file.exists(study_cache_file)) {
  study_info <- read_csv(study_cache_file, show_col_types = FALSE) %>% standardize_study_cols()
} else {
  study_info <- empty_study_tbl()
}

still_needed <- setdiff(unique_study_ids, unique(study_info$study_id))
message("  Studies already cached: ", length(unique_study_ids) - length(still_needed))
message("  Studies still needed: ",   length(still_needed))

if (length(still_needed) > 0) {
  for (i in seq_along(still_needed)) {
    sid      <- still_needed[i]
    new_info <- tryCatch(get_study_info(sid), error = function(e) NULL)
    if (!is.null(new_info)) {
      study_info <- bind_rows(study_info, standardize_study_cols(new_info)) %>%
        distinct(study_id, .keep_all = TRUE)
      write_csv(study_info, study_cache_file)
    }
    if (i %% 50 == 0 || i == length(still_needed))
      message("  Study info: ", i, "/", length(still_needed), " fetched")
  }
}

#========================
# Step 4: Join study metadata
#========================

gws <- gws %>% select(variant_id, association_id, study_id, pvalue, pvalue_mantissa, pvalue_exponent, queried_gene)

results <- if (nrow(study_info) > 0) {
  gws %>%
    left_join(study_info, by = "study_id") %>%
    arrange(queried_gene, pvalue) %>%
    distinct(queried_gene, study_id, reported_trait, .keep_all = TRUE)
} else {
  warning("No study info retrieved.")
  gws
}

#========================
# Step 5: Annotate reported traits with ChatGPT
#========================

message("Annotating traits with ChatGPT...")
trait_annotations <- annotate_traits_with_chatgpt_cached(results$reported_trait)
results <- results %>%
  select(-any_of(c("is_psychiatric", "disorder", "is_temperament", "temperament_category"))) %>%
  left_join(trait_annotations, by = c("reported_trait" = "trait"))

message("Annotation check:")
print(results %>% summarize(
  n_rows                = n(),
  n_missing_annotation  = sum(is.na(is_psychiatric)),
  n_psychiatric         = sum(is_psychiatric,  na.rm = TRUE),
  n_temperament         = sum(is_temperament,  na.rm = TRUE)
))

#========================
# Step 6: Output
#========================

message("Done. ", nrow(results), " GWS hits across ", n_distinct(results$queried_gene), " genes.")
message("  Psychiatric: ", sum(results$is_psychiatric, na.rm = TRUE))
message("  Temperament: ", sum(results$is_temperament, na.rm = TRUE))

results <- results %>% select(-any_of("temperament_category"))

write_csv(results, out_file)
message("Written to ", out_file)
