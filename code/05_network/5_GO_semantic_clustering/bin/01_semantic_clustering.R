#!/usr/bin/env Rscript
# ============================================================
# 01_semantic_clustering.R
#
# Build a semantically clustered, GPT-labelled table of GO enrichment that
# compares human-only ("run1" / ONLY) vs human+dogCD ("run2" / CCD) signal for
# DEP, OCD and SCH, plus the single-species, dog-only dogCD network (its own
# human_disease value, run1 only -- run2 is not applicable and left at 0),
# restricted to the top-level network community within each enrichment file.
#
# Pipeline:
#   1. Read per-disease GO enrichment files (run1 = human only, run2 = +dogCD),
#      plus the single dogCD enrichment file.
#   2. Within each file, pick the top-level community and keep only its terms.
#   3. Join run1/run2 per disease on the GO term; convert p-values to -log10.
#   4. Cluster significant GO terms by GOSemSim / Wang semantic similarity.
#   5. Score the clustering: per-cluster cohesion + representative term, and
#      a scan of alternative cutoffs (silhouette, cophenetic correlation).
#   6. Label each cluster with GPT (short label, wrapped label, description).
#   7. Write per-term and per-cluster tables (+ a Google Sheet tab).
#
# Inputs:
#   ../Figure4/<DISEASE>_{ONLY,CCD}_*_GO_enrichment.tsv  (see input_files)
#   ../Figure4/CCD_ONLY_hierachy_full_GO_enrichment.tsv  (dogCD, single run)
#
# Outputs (suffix encodes the clustering cutoff, e.g. "cut60"):
#   GO_clusters_labeled.<suffix>.txt                              (per cluster)
#   GO_with_without_dogCD.semantic_clustering.top_community.<suffix>.txt  (per term)
#   GO_cluster_cutoff_scan.top_community.txt                      (per cutoff)
#
# Caching: semantic clustering and GPT labels are cached so reruns resume
# rather than recompute / re-query the API.
# ============================================================

library(tidyverse)
library(glue)
library(httr2)
library(jsonlite)


#========================
# Package utilities
#========================

# GOSemSim / org.Hs.eg.db mask several tidyverse verbs, so they are loaded only
# for the clustering step and unloaded immediately afterwards. This helper
# detaches + unloads a package if (and only if) it is currently loaded.
unload_pkg_if_loaded <- function(pkg) {
  pkg_search_name <- paste0("package:", pkg)
  if (pkg_search_name %in% search()) {
    try(detach(pkg_search_name, unload = TRUE, character.only = TRUE), silent = TRUE)
  }
  if (pkg %in% loadedNamespaces()) {
    try(unloadNamespace(pkg), silent = TRUE)
  }
  invisible(NULL)
}

#========================
# Settings
#========================

GO_ontology    <- "BP"    # GO ontology used for semantic similarity (BP/MF/CC)
cluster_cutoff <- 0.6     # Wang similarity cutoff; clusters cut at height 1 - this
cutoff_scan    <- c(0.4, 0.5, 0.6, 0.7, 0.8)  # alternative cutoffs scored for comparison

gpt_model  <- "gpt-5"     # OpenAI model used for cluster labelling
batch_size <- 5           # clusters sent to GPT per API call
test_run   <- FALSE       # TRUE -> label only the first 2 clusters / 1 batch

# Network community column. Its name varies between files, so the first of
# these that is present in a given file is used.
community_cols <- c("community", "name.1")

# Keep only the single "top-level" community WITHIN EACH FILE. Community names
# differ between files, so the top community is chosen per file as the largest
# one: by query_size if that column exists, otherwise by number of GO-term rows.
# Set to FALSE to keep all communities (the community column is still retained).
restrict_to_top <- TRUE

# Filename suffix encoding the clustering cutoff (e.g. cluster_cutoff 0.6 -> "cut60").
suffix <- paste0("cut", round(cluster_cutoff * 100, 0))

# Cache files and final outputs. Both caches are durable: the clustering result
# (.rds) and the GPT cluster labels (.txt) persist across runs so neither the
# similarity matrix nor the GPT API is recomputed for work already done. Delete
# a cache file by hand to force that step to rerun.
cache_clustering <- paste0("cache_go_clusters.top_community.", suffix, ".rds")
cache_gpt        <- paste0("cache_gpt_labels.top_community.", suffix, ".txt")
# The similarity matrix does not depend on the cutoff, so it has no suffix.
cache_sim_mat    <- paste0("cache_go_sim_mat.top_community.", GO_ontology, ".rds")


out_file      <- paste0("GO_clusters_labeled.", suffix, ".txt")
out_file_long <- paste0("GO_with_without_dogCD.semantic_clustering.top_community.", suffix, ".txt")
out_file_scan <- "GO_cluster_cutoff_scan.top_community.txt"

#========================
# Exit early if final output already exists
#========================

# Avoid an expensive rerun (similarity matrix + GPT calls) if the per-term
# output and cutoff scan are already present. Delete one to force a fresh run.
if (file.exists(out_file_long) && file.exists(out_file_scan)) {
  message("Output already exists: ", out_file_long, ", ", out_file_scan)
  message("Skipping rerun.")
  stop("Done", call. = FALSE)
}

# One row per (disease, run) -> input enrichment file. run 1 = human only,
# run 2 = human + dogCD.
indir <- "./"  # input TSVs are staged into the working directory

input_files <- tibble(
  ds = c("DEP", "DEP", "OCD", "OCD", "SCH", "SCH"),
  run = c(1, 2, 1, 2, 1, 2),
  file = c(
    paste0(indir,"DEP_ONLY_hierarchy_full_GO_enrichment.tsv"),
    paste0(indir,"DEP_CCD_hierachy_full_GO_enrichment.tsv"),
    paste0(indir,"OCD_ONLY_hierachy_full_GO_enrichment.tsv"),
    paste0(indir,"Compulsive_hierachy_full_GO_enrichment_251023.tsv"),
    paste0(indir,"SCH_ONLY_hierachy_full_GO_enrichment.tsv"),
    paste0(indir,"SCH_CCD_hierachy_full_GO_enrichment.tsv")
  )
)

# Single-species, dog-only network: one enrichment file, no run1/run2 pair.
dogcd_file <- paste0(indir, "CCD_ONLY_hierachy_full_GO_enrichment.tsv")

#========================
# Read GO enrichment files
#========================

diseases <- c("DEP", "OCD", "SCH")

# Find the community column in a file (first match of community_cols).
find_community_col <- function(df, file) {
  hit <- intersect(community_cols, names(df))
  if (length(hit) == 0)
    stop("No community column (", paste(community_cols, collapse = ", "),
         ") found in ", file)
  hit[1]
}

# Pick the "top-level" community within one file: largest by query_size if that
# column is present, otherwise by the number of GO-term rows. As a sanity check,
# the largest community must ALSO be the lowest-numbered one (numeric suffix of
# the community name); otherwise we stop with an error.
top_community_of <- function(df, ccol, file = "") {
  comm <- df[[ccol]]
  if ("query_size" %in% names(df)) {
    # query_size is constant within a community; take its max per community.
    sz <- tapply(suppressWarnings(as.numeric(df$query_size)), comm,
                 function(v) max(v, na.rm = TRUE))
  } else {
    # fallback: the community contributing the most enriched GO-term rows.
    sz <- table(comm)
  }
  top_by_size <- names(sz)[which.max(sz)]
  
  # lowest-numbered community (trailing digits of the community name)
  comms <- names(sz)
  nums  <- suppressWarnings(as.numeric(stringr::str_extract(comms, "[0-9]+$")))
  if (all(is.na(nums)))
    stop("Could not extract community numbers from column '", ccol, "' in ", file)
  lowest_numbered <- comms[which.min(nums)]
  
  if (!identical(top_by_size, lowest_numbered))
    stop("Top community by size (", top_by_size, ") is not the lowest-numbered ",
         "community (", lowest_numbered, ") in ", file)
  
  top_by_size
}

# Read one enrichment file, optionally restrict to its top community, and return
# GO_id / GO_term / community / <p_out> (one row per term, smallest p kept).
# `p_out` names the p-value column in the result ("run1_p" or "run2_p").
read_enrichment <- function(file, p_out) {
  df   <- readr::read_tsv(file, show_col_types = FALSE)
  ccol <- find_community_col(df, file)
  if (restrict_to_top) {
    top <- top_community_of(df, ccol, file)
    message("  ", basename(file), ": top community = ", top,
            " (col '", ccol, "')")
    df <- df[df[[ccol]] == top, , drop = FALSE]
  }
  df %>%
    dplyr::transmute(GO_id = native, GO_term = name,
                     community = .data[[ccol]], p = p_value) %>%
    dplyr::group_by(GO_id, GO_term) %>%
    dplyr::summarize(!!p_out := min(p, na.rm = TRUE),
                     community = dplyr::first(community), .groups = "drop")
}

# For each disease, read its run1 (human only) and run2 (+dogCD) files and
# combine them into one table of -log10(p) signal per GO term.
raw_d <- purrr::map_dfr(diseases, function(in_disease) {
  
  only_file <- input_files %>%
    dplyr::filter(run == 1, ds == in_disease) %>%
    dplyr::pull(file)
  
  both_file <- input_files %>%
    dplyr::filter(run == 2, ds == in_disease) %>%
    dplyr::pull(file)
  
  if (length(only_file) != 1) stop("Missing or duplicated run1 file for ", in_disease)
  if (length(both_file) != 1) stop("Missing or duplicated run2 file for ", in_disease)
  
  # Community names differ between files, so the top community is chosen per
  # file and the run1/run2 tables are joined on the GO term only. The chosen
  # community is kept from each run (CCD -> community, human-only -> community_human).
  d1 <- read_enrichment(only_file, "run1_p") %>%
    dplyr::rename(community_human = community)
  d2 <- read_enrichment(both_file, "run2_p")
  
  # full_join keeps terms significant in either run; absent p-values become 1
  # (i.e. -log10(p) = 0) below.
  d2 %>%
    dplyr::full_join(d1, by = c("GO_id", "GO_term")) %>%
    dplyr::mutate(human_disease = in_disease) %>%
    replace_na(list(run2_p = 1, run1_p = 1))
})

# dogCD has a single enrichment file (no run1/run2 pair): run1 carries its
# -log10(p) signal, run2 is not applicable and set to 1 (-log10(p) = 0).
dogcd_d <- read_enrichment(dogcd_file, "run1_p") %>%
  dplyr::rename(community_human = community) %>%
  dplyr::mutate(community = community_human, human_disease = "dogCD", run2_p = 1)

raw_d <- dplyr::bind_rows(raw_d, dogcd_d)

# Convert p-values to -log10 signal and keep the columns used downstream.
df_annotated <- raw_d %>%
  dplyr::mutate(
    run1 = -log10(run1_p),   # human only
    run2 = -log10(run2_p)    # human + dogCD
  ) %>%
  dplyr::filter(!is.na(GO_id), !is.na(GO_term)) %>%
  dplyr::select(human_disease, community, community_human, GO_id, GO_term, run1, run2) %>%
  dplyr::arrange(human_disease, dplyr::desc(run2))

message("Retained GO term rows: ", nrow(df_annotated))

#========================
# Select terms to cluster
#========================

# Cluster only terms that are significant (-log10 p >= 5) in at least one run.
df_for_clustering <- df_annotated %>%
  dplyr::filter(run1 >= 5 | run2 >= 5) %>%
  dplyr::select(GO_id, GO_term) %>%
  dplyr::distinct()

message("Unique GO terms selected for clustering: ", nrow(df_for_clustering))

#========================
# Semantic clustering
#========================

# Cluster GO terms by Wang semantic similarity (average-linkage hclust cut at
# height 1 - cluster_cutoff). Result cached to avoid recomputing the matrix.
if (file.exists(cache_clustering)) {
  
  message("Loading clustering result from cache: ", cache_clustering)
  df_with_clusters <- readRDS(cache_clustering)
  
} else {
  
  message("Loading GOSemSim only for semantic clustering.")
  requireNamespace("GOSemSim", quietly = TRUE)
  requireNamespace("org.Hs.eg.db", quietly = TRUE)
  
  message("Loading GOSemSim semantic data for ontology: ", GO_ontology)
  
  # Precompute the GO DAG / IC-free semantic data for the chosen ontology.
  sem_data <- GOSemSim::godata(
    annoDb = "org.Hs.eg.db",
    ont = GO_ontology,
    computeIC = FALSE
  )
  
  dat <- df_for_clustering
  ids <- unique(dat$GO_id)
  
  if (length(ids) == 0) {
    stop("No GO terms passed the clustering filter.")
  }
  
  message("Computing semantic similarity matrix for ", length(ids), " GO terms...")
  
  # Pairwise Wang similarity (NA -> 0 so dissimilarity is well defined).
  sim_mat <- GOSemSim::mgoSim(
    ids,
    ids,
    sem_data,
    measure = "Wang",
    combine = NULL
  )
  
  sim_mat[is.na(sim_mat)] <- 0
  saveRDS(sim_mat, cache_sim_mat)
  
  if (length(ids) == 1) {
    # Degenerate case: a single term is its own cluster.
    cluster_df <- tibble(
      GO_id = ids,
      GO_semantic_cluster = "cluster_1"
    )
  } else {
    # Average-linkage hierarchical clustering on (1 - similarity), cut at the
    # cutoff to assign each GO term to a semantic cluster.
    hc <- stats::hclust(
      stats::as.dist(1 - sim_mat),
      method = "average"
    )
    
    clusters <- stats::cutree(
      hc,
      h = 1 - cluster_cutoff
    )
    
    cluster_df <- tibble(
      GO_id = names(clusters),
      GO_semantic_cluster = paste0("cluster_", clusters)
    )
  }
  
  # Attach cluster IDs back to every annotated row (terms not clustered, i.e.
  # below the significance filter, get NA).
  df_clustered_terms <- dat %>%
    dplyr::left_join(cluster_df, by = "GO_id") %>%
    dplyr::distinct()
  
  df_with_clusters <- df_annotated %>%
    dplyr::left_join(
      df_clustered_terms %>%
        dplyr::select(GO_id, GO_semantic_cluster) %>%
        dplyr::distinct(),
      by = "GO_id"
    )
  
  saveRDS(df_with_clusters, cache_clustering)
  message("Saved clustering cache: ", cache_clustering)
  
  unload_pkg_if_loaded("GOSemSim")
  unload_pkg_if_loaded("org.Hs.eg.db")
  message("Unloaded GOSemSim/org.Hs.eg.db after semantic clustering.")
}

#========================
# Cluster quality metrics
#========================

# Quantitative support for the clustering, from the same Wang similarity
# matrix used to build it (cached; recomputed only if missing or stale).
metric_ids <- unique(df_for_clustering$GO_id)

sim_mat <- if (file.exists(cache_sim_mat)) readRDS(cache_sim_mat) else NULL

if (is.null(sim_mat) || !all(metric_ids %in% rownames(sim_mat))) {
  message("Computing semantic similarity matrix for cluster metrics...")
  requireNamespace("GOSemSim", quietly = TRUE)
  requireNamespace("org.Hs.eg.db", quietly = TRUE)
  sem_data <- GOSemSim::godata(annoDb = "org.Hs.eg.db", ont = GO_ontology,
                               computeIC = FALSE)
  sim_mat <- GOSemSim::mgoSim(metric_ids, metric_ids, sem_data,
                              measure = "Wang", combine = NULL)
  sim_mat[is.na(sim_mat)] <- 0
  saveRDS(sim_mat, cache_sim_mat)
  unload_pkg_if_loaded("GOSemSim")
  unload_pkg_if_loaded("org.Hs.eg.db")
}
sim_mat <- sim_mat[metric_ids, metric_ids, drop = FALSE]

# Rebuild the tree exactly as in the clustering step.
sim_dist <- stats::as.dist(1 - sim_mat)
hc       <- stats::hclust(sim_dist, method = "average")

# Sanity check: cutting the rebuilt tree at cluster_cutoff must reproduce the
# cached cluster assignments (same partition; cluster numbers may differ).
cluster_assign <- df_with_clusters %>%
  dplyr::filter(!is.na(GO_semantic_cluster)) %>%
  dplyr::distinct(GO_id, .keep_all = TRUE) %>%
  dplyr::select(GO_id, GO_term, GO_semantic_cluster)

recut   <- stats::cutree(hc, h = 1 - cluster_cutoff)[cluster_assign$GO_id]
xtab    <- table(cluster_assign$GO_semantic_cluster, recut) > 0
if (!all(rowSums(xtab) == 1) || !all(colSums(xtab) == 1)) {
  stop("Rebuilt tree does not reproduce the cached clusters; delete ",
       cache_clustering, " and ", cache_sim_mat, " and rerun.")
}

# Per cluster: cohesion (mean and minimum pairwise Wang similarity among its
# terms; NA for single-term clusters) and a representative term, the one with
# the highest mean similarity to the rest of its cluster.
cluster_metrics <- cluster_assign %>%
  dplyr::group_by(GO_semantic_cluster) %>%
  dplyr::group_modify(function(g, key) {
    m      <- sim_mat[g$GO_id, g$GO_id, drop = FALSE]
    pairs  <- m[upper.tri(m)]
    rep_i  <- which.max(rowMeans(m))
    tibble(
      within_cluster_mean_sim = if (length(pairs)) round(mean(pairs), 3) else NA_real_,
      within_cluster_min_sim  = if (length(pairs)) round(min(pairs), 3) else NA_real_,
      representative_GO_id    = g$GO_id[rep_i],
      representative_GO_term  = g$GO_term[rep_i]
    )
  }) %>%
  dplyr::ungroup()

# Cutoff scan: how the clustering changes across candidate cutoffs. Mean
# silhouette width (-1 to 1) measures how much closer terms are to their own
# cluster than to the nearest other cluster; single-term clusters score 0.
# Cophenetic correlation measures how well the tree preserves the original
# similarities (one value for the tree, so the same at every cutoff).
cophenetic_cor <- stats::cor(sim_dist, stats::cophenetic(hc))

cutoff_scan_df <- purrr::map_dfr(sort(unique(c(cutoff_scan, cluster_cutoff))), function(co) {
  cl    <- stats::cutree(hc, h = 1 - co)
  sizes <- table(cl)
  sil   <- if (length(sizes) > 1 && length(sizes) < length(cl)) {
    mean(cluster::silhouette(cl, sim_dist)[, "sil_width"])
  } else NA_real_
  # Cohesion averaged over multi-term clusters.
  within <- purrr::map_dbl(names(sizes)[sizes > 1], function(k) {
    m <- sim_mat[cl == as.integer(k), cl == as.integer(k)]
    mean(m[upper.tri(m)])
  })
  tibble(
    cutoff                  = co,
    used                    = isTRUE(all.equal(co, cluster_cutoff)),
    n_terms                 = length(cl),
    n_clusters              = length(sizes),
    n_singletons            = sum(sizes == 1),
    median_cluster_size     = stats::median(sizes),
    max_cluster_size        = max(sizes),
    mean_within_cluster_sim = round(mean(within), 3),
    mean_silhouette         = round(sil, 3),
    cophenetic_cor          = round(cophenetic_cor, 3)
  )
})

readr::write_tsv(cutoff_scan_df, out_file_scan)
message("Wrote: ", out_file_scan)
print(cutoff_scan_df)

#========================
# Build cluster-term membership
#========================

# Distinct (cluster, GO term) memberships, then one row per cluster with its
# concatenated terms/ids and a stable label_id for the GPT batches.
cluster_terms <- df_with_clusters %>%
  dplyr::filter(!is.na(GO_semantic_cluster), !is.na(GO_term)) %>%
  dplyr::select(GO_semantic_cluster, GO_id, GO_term) %>%
  dplyr::distinct()

clusters_to_label <- cluster_terms %>%
  dplyr::group_by(GO_semantic_cluster) %>%
  dplyr::summarize(
    GO_terms = paste(unique(GO_term), collapse = "; "),
    GO_ids = paste(unique(GO_id), collapse = "; "),
    n_terms = dplyr::n_distinct(GO_term),
    .groups = "drop"
  ) %>%
  dplyr::arrange(GO_semantic_cluster) %>%
  dplyr::mutate(label_id = dplyr::row_number())

if (test_run) {
  message("test_run=TRUE: limiting to 2 clusters")
  clusters_to_label <- clusters_to_label %>%
    dplyr::slice_head(n = 2)
}

message("Total clusters to label: ", nrow(clusters_to_label))

#========================
# GPT labeling
#========================

if (Sys.getenv("OPENAI_API_KEY") == "") {
  stop("OPENAI_API_KEY is not set. Run: Sys.setenv(OPENAI_API_KEY = 'your_key_here')")
}

`%||%` <- function(x, y) if (is.null(x)) y else x

# Flatten the OpenAI Responses API output into a single text string.
extract_response_text <- function(out) {
  text_blocks <- unlist(
    lapply(out$output, function(x) {
      if (!is.null(x$content)) {
        sapply(x$content, function(y) y$text %||% "")
      } else {
        NULL
      }
    })
  )
  paste(text_blocks, collapse = "\n")
}

# Strip ```json ... ``` fences GPT sometimes wraps around JSON.
strip_json_fences <- function(x) {
  x <- stringr::str_trim(x)
  x <- stringr::str_remove(x, "^```json\\s*")
  x <- stringr::str_remove(x, "^```\\s*")
  x <- stringr::str_remove(x, "\\s*```$")
  stringr::str_trim(x)
}

# Call the OpenAI Responses API with retries + backoff; returns the text output.
run_gpt_prompt <- function(system_prompt, user_prompt, model = "gpt-5", max_tries = 3) {
  
  for (try_i in seq_len(max_tries)) {
    
    message("    API attempt ", try_i, "/", max_tries)
    
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
      httr2::req_error(is_error = function(resp) FALSE) %>%   # handle status manually
      httr2::req_perform()
    
    if (httr2::resp_status(resp) < 400) {
      return(extract_response_text(httr2::resp_body_json(resp)))
    }
    
    message("    API error: ", httr2::resp_body_string(resp))
    
    if (try_i < max_tries) {
      Sys.sleep(5 * try_i)          # linear backoff between attempts
    } else {
      stop(httr2::resp_body_string(resp))
    }
  }
}

message("\n--- GPT labeling ---")

# Resume from a checkpoint if present: skip clusters already labelled.
if (file.exists(cache_gpt)) {
  
  message("Resuming from GPT checkpoint: ", cache_gpt)
  
  existing_labels <- readr::read_tsv(cache_gpt, show_col_types = FALSE)
  
  clusters_to_label <- clusters_to_label %>%
    dplyr::anti_join(
      existing_labels %>% dplyr::select(GO_semantic_cluster),
      by = "GO_semantic_cluster"
    )
  
  message("Remaining clusters after checkpoint: ", nrow(clusters_to_label))
  
} else {
  
  existing_labels <- tibble(
    GO_semantic_cluster = character(),
    cluster_label = character(),
    cluster_label_wrapped = character(),
    cluster_description = character()
  )
}

# Instructions defining the JSON contract GPT must return for each cluster.
system_prompt <- "
You are an expert computational biologist.

For each group of semantically related Gene Ontology terms, generate:
1. a short cluster label of 4-8 words
2. a display-formatted cluster label constrained to a maximum line width of 20 characters
3. a one-sentence biological description

Be concise and Nature Results-style.

For cluster_label:
- Use plain text with no formatting.

For cluster_label_wrapped:
- Format the label so no line exceeds 18-22 characters.
- Represent line breaks using exactly two underscores: __
- Add hyphens only when needed to split long words across lines.
- Use linguistically sensible hyphenation (hyphenate between syllables).
- Use biologically sensible hyphenation.
- Replace `and` with `&`
- Keep wording concise and readable.
- Do not use actual newline characters.

Example:
Synaptic vesicle trafficking
could become:
Synaptic vesicle__trafficking

Another example:
Synaptic modulation and plasticity
could become:
Synaptic modul-__ation and plasticity

Return valid JSON only.
Return an array of objects with exactly these fields:
label_id, cluster_label, cluster_label_wrapped, cluster_description
"

# Label clusters in batches, checkpointing after each batch so a failure or
# interruption can resume without re-querying already-labelled clusters.
if (nrow(clusters_to_label) > 0) {
  
  batches <- split(
    clusters_to_label,
    ceiling(seq_len(nrow(clusters_to_label)) / batch_size)
  )
  
  for (batch_i in seq_along(batches)) {
    
    if (test_run && batch_i > 1) {
      break
    }
    
    batch <- batches[[batch_i]]
    
    message(
      "\nBatch ", batch_i, "/", length(batches),
      " | label_id ", min(batch$label_id), "-", max(batch$label_id),
      " | n=", nrow(batch)
    )
    
    # Send only the fields GPT needs; it returns labels keyed by label_id.
    batch_payload <- batch %>%
      dplyr::select(label_id, GO_semantic_cluster, n_terms, GO_terms)
    
    user_prompt <- glue::glue(
      "Label the following GO semantic clusters.\n\n{jsonlite::toJSON(batch_payload, dataframe = 'rows', auto_unbox = TRUE, pretty = TRUE)}"
    )
    
    gpt_text <- run_gpt_prompt(
      system_prompt = system_prompt,
      user_prompt = user_prompt,
      model = gpt_model
    )
    
    parsed <- tryCatch(
      jsonlite::fromJSON(strip_json_fences(gpt_text)) %>%
        as_tibble(),
      error = function(e) {
        message("Could not parse GPT output as JSON:\n", gpt_text)
        stop(e)
      }
    )
    
    # Re-key GPT's labels back onto this batch's clusters via label_id.
    batch_labels <- batch %>%
      dplyr::select(label_id, GO_semantic_cluster) %>%
      dplyr::left_join(parsed, by = "label_id") %>%
      dplyr::mutate(
        cluster_label = trimws(cluster_label),
        cluster_label_wrapped = trimws(cluster_label_wrapped),
        cluster_description = trimws(cluster_description)
      ) %>%
      dplyr::select(
        GO_semantic_cluster,
        cluster_label,
        cluster_label_wrapped,
        cluster_description
      )
    
    existing_labels <- dplyr::bind_rows(existing_labels, batch_labels) %>%
      dplyr::distinct(GO_semantic_cluster, .keep_all = TRUE)
    
    readr::write_tsv(existing_labels, cache_gpt)
    message("  Checkpoint saved: ", cache_gpt)
  }
  
} else {
  message("No new clusters to label.")
}

#========================
# Write final output
#========================

message("\n--- Writing output ---")

# Prefer the on-disk checkpoint (authoritative, de-duplicated) over the
# in-memory labels.
if (file.exists(cache_gpt)) {
  cluster_labels <- readr::read_tsv(cache_gpt, show_col_types = FALSE) %>%
    dplyr::distinct(GO_semantic_cluster, .keep_all = TRUE) %>%
    dplyr::arrange(GO_semantic_cluster)
} else {
  cluster_labels <- existing_labels
}

# Per-cluster summary: terms + cohesion metrics + GPT labels (one row per cluster).
cluster_summary_labeled <- cluster_terms %>%
  dplyr::group_by(GO_semantic_cluster) %>%
  dplyr::summarize(
    GO_terms = paste(unique(GO_term), collapse = "; "),
    GO_ids = paste(unique(GO_id), collapse = "; "),
    n_terms = dplyr::n_distinct(GO_term),
    .groups = "drop"
  ) %>%
  dplyr::left_join(cluster_metrics, by = "GO_semantic_cluster") %>%
  dplyr::left_join(cluster_labels, by = "GO_semantic_cluster") %>%
  dplyr::arrange(GO_semantic_cluster)

readr::write_tsv(cluster_summary_labeled, out_file)
message("Wrote: ", out_file)

# Per-term table: every GO term with its run1/run2 signal, communities, cluster,
# and the cluster's GPT labels.
df_out <- df_with_clusters %>%
  dplyr::left_join(cluster_labels, by = "GO_semantic_cluster")

readr::write_tsv(df_out, out_file_long)
message("Wrote: ", out_file_long)

# NOTE: cache_gpt is kept (NOT deleted) so it acts as a durable label cache.
# On any rerun, clusters already present here are skipped via the anti_join in
# the GPT-labeling step, so GPT is only queried for genuinely new clusters.
# Delete this file by hand to force every cluster to be re-labelled.
message("Kept GPT label cache: ", cache_gpt)

unload_pkg_if_loaded("GOSemSim")
unload_pkg_if_loaded("org.Hs.eg.db")
message("\nDone.")