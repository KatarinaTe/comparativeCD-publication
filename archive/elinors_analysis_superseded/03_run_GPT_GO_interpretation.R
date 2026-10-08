# ============================================================
# 03_run_GPT_GO_interpretation.R
#
# Purpose:
#   Run ChatGPT interpretation of GO enrichment results using
#   only semantic-clustered GO results.
#
# This script assumes that 01_prepare_GO_tables.R creates:
#   idf_for_gpt_go
#   summary_for_gpt
#   cluster_summary_for_gpt
#
# Required columns include:
#   GO_id
#   GO_ontology
#   GO_semantic_cluster
#
# GO_semantic_cluster should be generated using GOSemSim-based
# ontology-aware semantic similarity clustering.
# ============================================================

library(glue)
library(jsonlite)
library(httr2)
library(readr)
library(tibble)
library(stringr)

#========================
# Load prepared GO tables
#========================

source("01_prepare_GO_tables.R")
#source("02_label_GO_clusters_with_GPT.R")
#========================
# Check required objects
#========================

required_objects <- c(
  "df_for_gpt_go",
  "summary_for_gpt",
  "cluster_summary_for_gpt"
)

missing_objects <- required_objects[!vapply(required_objects, exists, logical(1))]

if (length(missing_objects) > 0) {
  stop(
    "Missing required object(s) from 01_prepare_GO_tables.R: ",
    paste(missing_objects, collapse = ", ")
  )
}

required_cols <- c(
  "GO_id",
  "GO_ontology",
  "GO_semantic_cluster"
)

missing_cols <- setdiff(required_cols, colnames(df_for_gpt_go))

if (length(missing_cols) > 0) {
  stop(
    "df_for_gpt_go is missing required column(s): ",
    paste(missing_cols, collapse = ", ")
  )
}
#========================
# Keep GPT input focused on significant results
#========================
# Retain clusters/GO terms if they exceed significance
# in at least one disease.

cluster_summary_for_gpt <- cluster_summary_for_gpt %>%
  dplyr::group_by(GO_semantic_cluster) %>%
  dplyr::filter(any(max_run1 > 5 | max_run2 > 5, na.rm = TRUE)) %>%
  dplyr::ungroup()

df_for_gpt_go <- df_for_gpt_go %>%
  dplyr::group_by(GO_semantic_cluster) %>%
  dplyr::filter(any(run1 > 5 | run2 > 5, na.rm = TRUE)) %>%
  dplyr::ungroup()

message("GPT cluster rows: ", nrow(cluster_summary_for_gpt))
message("GPT GO-term rows: ", nrow(df_for_gpt_go))

#========================
# Helper functions
#========================

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

extract_response_text <- function(out) {
  text_out <- unlist(
    lapply(out$output, function(x) {
      if (!is.null(x$content)) {
        sapply(x$content, function(y) y$text %||% "")
      } else {
        NULL
      }
    })
  )
  
  paste(text_out, collapse = "\n")
}

#========================
# Check API key
#========================

if (Sys.getenv("OPENAI_API_KEY") == "") {
  stop(
    "OPENAI_API_KEY is not set. Run: ",
    "Sys.setenv(OPENAI_API_KEY = 'your_key_here')"
  )
}

#========================
# OpenAI API wrapper
#========================

run_gpt_prompt <- function(system_prompt, user_prompt, model = "gpt-5") {
  
  resp <- httr2::request("https://api.openai.com/v1/responses") %>%
    httr2::req_headers(
      Authorization = paste("Bearer", Sys.getenv("OPENAI_API_KEY")),
      `Content-Type` = "application/json"
    ) %>%
    httr2::req_body_json(
      list(
        model = model,
        input = list(
          list(role = "system", content = system_prompt),
          list(role = "user", content = user_prompt)
        )
      )
    ) %>%
    httr2::req_error(is_error = function(resp) FALSE) %>%
    httr2::req_perform()
  
  if (httr2::resp_status(resp) >= 400) {
    stop(httr2::resp_body_string(resp))
  }
  
  out <- httr2::resp_body_json(resp)
  
  extract_response_text(out)
}

#========================
# Semantic-clustered interpretation
#========================

system_prompt_cluster <- "
You are an expert computational biologist and statistical genomics analyst.

Interpret Gene Ontology enrichment results for psychiatric disease gene networks.

run1 is the human GWAS-only gene network.
run2 is the colocalized gene network from human GWAS and gene network from dog GWAS from NetColoc.

Scores are -log10(p-values); values >5 are significant.

Use GO_id, GO_ontology, and GO_semantic_cluster to account for GO hierarchy and redundancy.

GO_semantic_cluster groups related GO terms using ontology-aware semantic similarity. Terms in the same cluster should generally be interpreted as one biological module rather than as independent findings.

Focus on biological themes and differences between run1 and run2.

Write a human-readable summary.

Use disease names rather than acronyms:
DEP = major depression
SCH = schizophrenia
OCD = OCD
"

user_prompt_cluster <- glue::glue("
Analyze the GO enrichment results below.

Definitions:
- run1 = gene network from human GWAS only
- run2 = colocalized gene network from human GWAS and gene network from dog GWAS
- Scores are -log10(p-values)
- A score >5 is significant
- delta = run2 - run1
- GO_id and GO_ontology provide Gene Ontology hierarchy context
- GO_semantic_cluster is a semantic redundancy cluster generated using GOSemSim/Wang similarity

Context:
- Major depression and schizophrenia are based on large, well-powered GWAS.
- OCD GWAS is smaller and less powered.
- Determine how adding dog GWAS impacts results for each disease

Return:
- Compare run2 with run1 across diseases.
- Quantify newly confident terms.
- Use GO_semantic_cluster to collapse redundant GO terms into biological modules.
- Do not interpret individual GO terms as independent evidence when they fall within the same semantic cluster.
- Describe any clusters that are more significant for OCD than other diseases
- Describe clusters that lose significances in colocalized network.
- A coherent interpretation.

Summary table:
{jsonlite::toJSON(summary_for_gpt, dataframe = 'rows', pretty = FALSE, auto_unbox = TRUE)}

Semantic cluster summary:
{jsonlite::toJSON(cluster_summary_for_gpt, dataframe = 'rows', pretty = FALSE, auto_unbox = TRUE)}

GO term table:
{jsonlite::toJSON(df_for_gpt_go, dataframe = 'rows', pretty = FALSE, auto_unbox = TRUE)}
")

interpret_cluster <- run_gpt_prompt(
  system_prompt = system_prompt_cluster,
  user_prompt = user_prompt_cluster,
  model = "gpt-5"
)

#========================
# Save output
#========================

out_file <- "GO_interpret_network_overlap_semantic_clustering.md"

readr::write_lines(
  c(
    "# GO interpretation using semantic clustering",
    "",
    interpret_cluster
  ),
  out_file
)

message("Wrote: ", normalizePath(out_file))

browseURL(normalizePath(out_file))

out_file <- "GO_interpret_network_overlap_semantic_clustering.txt"

readr::write_lines(
  c(
    "# GO interpretation using semantic clustering",
    "",
    interpret_cluster
  ),
  out_file
)

message("Wrote: ", normalizePath(out_file))
