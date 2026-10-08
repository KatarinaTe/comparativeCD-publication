#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# plot_gene_enrichment_grids.v4.R
#
# Changes vs v3:
#   1. p_cutoff = 1e-5 (was 1e-4).
#   2. Title shown only on the FIRST panel of each community; subsequent
#      panels have no title, and their height is reduced by title_pad_in.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(ape)
  library(cowplot)
})
stopifnot(requireNamespace("GO.db",        quietly = TRUE),
          requireNamespace("AnnotationDbi", quietly = TRUE))

`%||%` <- function(a, b) if (is.null(a)) b else a

## ---------------------------------------------------------------------------
## CONFIG
## ---------------------------------------------------------------------------
cache_path <- "all_go_combined.rds"
plot_ds    <- "OCD_CCD"

# Communities: same five as GO_enrichment_OCD_CCD_facets.p5.barchart.pdf
communities <- c("C185", "C187", "C186", "C197")

# --- barchart / focal-term parameters (match plot_fold_enrichment.R) --------
# These define which GO terms are focal (★, bolded). They must match the
# filters used by plot_fold_enrichment.R to produce the barchart PDF.
barchart_pvalue_max    <- 1e-5
barchart_max_sets      <- 10
barchart_min_intersect <- 5

# --- pool of terms eligible to appear as rows (focals + their relatives) ----
# Only terms with p < p_cutoff can appear as ancestors / siblings of focals.
p_cutoff        <- 1e-5
max_termsize    <- 1000
min_term_size   <- 5
min_overlap     <- 3

include_ancestors <- TRUE
include_siblings  <- TRUE
include_children  <- FALSE

# --- merging focal terms into fewer panels per community --------------------
merge_focals    <- TRUE
merge_dist      <- 3
max_panel_terms <- 25

# --- seed genes (frozen CSV) ------------------------------------------------
seed_file   <- "gene_list_gwas_catalog_overlap_with_region.csv"  # frozen export of the seed-gene sheet tab below
seed_tab    <- "gene list w gwas catalog overlap"
col_gene_id <- "gene"
col_dogseed <- "dog seed"
col_humseed <- "human seed"
col_dogcdp  <- "dogCD P"          # p-value column for dogCD GWAS
dogcd_p_cutoff <- 4e-7            # significance threshold

# --- authoritative GO annotations -------------------------------------------
goa_gaf <- "goa_human.gaf.gz"
goa_url <- "http://current.geneontology.org/annotations/goa_human.gaf.gz"

# --- layout -----------------------------------------------------------------
square_tiles  <- TRUE
ylab_char_in  <- 0.048
row_height_in <- 0.15
# Panel height formula: (n_rows + y_data_pad) * row_height_in * height_scale + axis_oh_in
# y_data_pad: extra y-range in data units beyond the data rows, from bar-axis ticks
#   (tick labels at min(tip_y)-1.05) plus expand(add=c(1.5,0.5)) → 2.55 units.
# axis_oh_in: fixed ggplot overhead from the top axis-text labels (gene names, inches).
# Together these ensure every panel has the same physical row height.
y_data_pad    <- 2.55
axis_oh_in    <- 0.30
height_scale  <- 0.75
tree_scale    <- 0.8
tree_gap      <- 0.2
bar_w_units   <- 6
bar_gap       <- 0.8

# --- column names -----------------------------------------------------------
col_ds       <- NULL
col_comm     <- "community"
col_go       <- "native"
col_name     <- "name"
col_sig      <- "significant"
col_pval     <- "p_value"
col_termsize <- "term_size"
col_qsize    <- "query_size"
col_source   <- "source"
col_genes    <- NULL
col_effdom   <- NULL

# --- community label colours (one per community) ----------------------------
comm_colors <- c(C185 = "#4e79a7",   # parent
                 C187 = "#9bbdd4",   # subcommunity of C185 (lighter blue)
                 C186 = "#59a14f",
                 C197 = "#e15759")

# --- gene-category colours --------------------------------------------------
cat_levels <- c("network gene", "dog seed", "human seed",
                "network & dog seed", "network & human seed")
cat_colors <- c("network gene"         = "#c39bd3",
                "dog seed"             = "#e1706c",
                "human seed"           = "#5c80bd",
                "network & dog seed"   = "#9e3b34",
                "network & human seed" = "#36568c")

## ---------------------------------------------------------------------------
## LOAD SEED GENES
## ---------------------------------------------------------------------------

as_lgl_flag <- function(x) {
  if (is.logical(x)) return(x %in% TRUE)
  if (is.numeric(x)) return(!is.na(x) & x != 0)
  toupper(trimws(as.character(x))) %in% c("TRUE", "T", "YES", "Y", "1", "X")
}

message("Reading seed genes from ", seed_file, " (tab '", seed_tab, "') ...")
seed_tbl <- utils::read.csv(seed_file, check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(all(c(col_gene_id, col_dogseed, col_humseed) %in% names(seed_tbl)))

norm_sym <- function(x) toupper(trimws(as.character(x)))
dog_seed_genes   <- unique(norm_sym(seed_tbl[[col_gene_id]][as_lgl_flag(seed_tbl[[col_dogseed]])]))
human_seed_genes <- unique(norm_sym(seed_tbl[[col_gene_id]][as_lgl_flag(seed_tbl[[col_humseed]])]))
dog_seed_genes   <- dog_seed_genes[!is.na(dog_seed_genes)   & dog_seed_genes   != "NA" & dog_seed_genes   != ""]
human_seed_genes <- human_seed_genes[!is.na(human_seed_genes) & human_seed_genes != "NA" & human_seed_genes != ""]
both_seed <- intersect(dog_seed_genes, human_seed_genes)
if (length(both_seed))
  message("Note: ", length(both_seed), " gene(s) flagged as BOTH dog and human seed.")
message("Seed genes: ", length(dog_seed_genes), " dog, ", length(human_seed_genes), " human")

# Dog seed genes with significant dogCD GWAS association (p ≤ dogcd_p_cutoff).
# These tiles get an asterisk overlay.
if (col_dogcdp %in% names(seed_tbl)) {
  pvals <- suppressWarnings(as.numeric(seed_tbl[[col_dogcdp]]))
  dogcd_sig_genes <- unique(norm_sym(
    seed_tbl[[col_gene_id]][!is.na(pvals) & pvals <= dogcd_p_cutoff]
  ))
  dogcd_sig_genes <- dogcd_sig_genes[dogcd_sig_genes != "" & dogcd_sig_genes != "NA"]
  message("dogCD significant genes (p <= ", dogcd_p_cutoff, "): ",
          length(dogcd_sig_genes), " — ", paste(dogcd_sig_genes, collapse = ", "))
} else {
  warning("Column '", col_dogcdp, "' not found in seed sheet — dogCD asterisks disabled.")
  dogcd_sig_genes <- character(0)
}

## ---------------------------------------------------------------------------
## GO-TERM -> human gene SYMBOLS (authoritative GOA annotations, DIRECT only)
## ---------------------------------------------------------------------------
if (!file.exists(goa_gaf)) {
  message("Downloading authoritative GO annotations:\n  ", goa_url)
  utils::download.file(goa_url, goa_gaf, mode = "wb", quiet = TRUE)
}
message("Reading GO annotations from ", goa_gaf, " ...")
gaf <- readr::read_tsv(goa_gaf, comment = "!", col_names = FALSE, quote = "",
                       col_types = readr::cols(.default = "c"), progress = FALSE)
gaf_direct <- gaf %>%
  dplyr::filter(X9 %in% c("P", "F", "C"),
                !grepl("NOT", X4, fixed = TRUE)) %>%
  dplyr::transmute(symbol = norm_sym(X3), go = X5) %>%
  dplyr::distinct()
direct2genes <- split(gaf_direct$symbol, gaf_direct$go)
message("GOA: ", nrow(gaf_direct), " direct annotations across ",
        length(direct2genes), " GO terms")

genes_in_term <- local({
  cache <- new.env(parent = emptyenv())
  function(go_id) {
    if (is.null(go_id) || is.na(go_id) || !grepl("^GO:", go_id))
      return(character(0))
    if (!is.null(cache[[go_id]])) return(cache[[go_id]])
    syms <- unique(direct2genes[[go_id]] %||% character(0))
    syms <- syms[!is.na(syms)]
    cache[[go_id]] <- syms
    syms
  }
})

## ---------------------------------------------------------------------------
## LOAD + INSPECT
## ---------------------------------------------------------------------------
stopifnot(file.exists(cache_path))
go <- readRDS(cache_path)
message("Loaded ", nrow(go), " rows x ", ncol(go), " cols")
glimpse(go)

## ---------------------------------------------------------------------------
## COLUMN DETECTION
## ---------------------------------------------------------------------------
detect_ds_col <- function(df, val) {
  cand <- names(df)[map_lgl(df, is.character)]
  hit  <- cand[map_lgl(cand, ~ any(df[[.x]] == val, na.rm = TRUE))]
  if (length(hit) == 0)
    stop("No character column contains '", val, "'. Set col_ds manually.")
  hit[1]
}

detect_genes_col <- function(df, exclude = c("evidences", "parents")) {
  cand <- setdiff(names(df)[map_lgl(df, is.character)], exclude)
  looks_list <- map_lgl(cand, function(nm) {
    x <- df[[nm]][!is.na(df[[nm]])]
    length(x) > 0 && mean(str_detect(head(x, 50), "^\\s*\\[")) > 0.5
  })
  hit <- cand[looks_list]
  if (length(hit) == 0)
    stop("Could not auto-detect the gene-list column. Set col_genes manually.")
  if (length(hit) > 1) {
    score <- map_dbl(hit, function(nm) {
      toks <- unlist(str_extract_all(head(df[[nm]], 20), "[A-Za-z0-9._-]+"))
      if (length(toks) == 0) return(0)
      mean(str_detect(toks, "^[A-Z0-9.-]+$"))
    })
    hit <- hit[which.max(score)]
  }
  hit[1]
}

detect_effdom_col <- function(df, exclude) {
  num <- setdiff(names(df)[map_lgl(df, is.numeric)], exclude)
  ok  <- num[map_lgl(num, function(nm) {
    v <- df[[nm]][!is.na(df[[nm]])]
    length(v) > 0 && min(v) > 1000 && dplyr::n_distinct(v) <= 12
  })]
  if (length(ok) == 0) return(NA_character_)
  ok[which.max(map_dbl(ok, ~ stats::median(df[[.x]], na.rm = TRUE)))]
}

col_ds    <- col_ds    %||% detect_ds_col(go, plot_ds)
col_genes <- col_genes %||% detect_genes_col(go)
message("ds column: ", col_ds, "  |  genes column: ", col_genes)

## ---------------------------------------------------------------------------
## FILTER + DERIVE
## ---------------------------------------------------------------------------
parse_genes <- function(x) {
  x %>%
    str_remove_all("[\\[\\]'\"]") %>%
    str_split(",\\s*") %>%
    map(~ str_trim(.x)) %>%
    map(~ .x[.x != "" & !is.na(.x)])
}

go_ds <- go %>%
  dplyr::filter(.data[[col_ds]]       == plot_ds,
                .data[[col_sig]]      == TRUE,
                .data[[col_pval]]     <  p_cutoff,
                .data[[col_termsize]] >  min_term_size,
                .data[[col_termsize]] <= max_termsize)
message(nrow(go_ds), " eligible rows (significant, p < ", p_cutoff,
        ", ", min_term_size, " < term_size <= ", max_termsize, ") for ds=", plot_ds)
stopifnot(nrow(go_ds) > 0)

col_effdom <- col_effdom %||%
  detect_effdom_col(go_ds, exclude = c(col_termsize, col_qsize))
if (is.na(col_effdom))
  stop("Could not find an effective_domain_size column; set col_effdom manually.")
message("effective_domain_size column: ", col_effdom)

go_ds <- go_ds %>%
  dplyr::mutate(
    .comm   = .data[[col_comm]],
    .go     = .data[[col_go]],
    .name   = .data[[col_name]],
    .source = .data[[col_source]],
    .pval   = .data[[col_pval]],
    .genes  = parse_genes(.data[[col_genes]]),
    .isize  = lengths(.genes),
    .fe     = (.isize / .data[[col_qsize]]) /
      (.data[[col_termsize]] / .data[[col_effdom]])
  ) %>%
  dplyr::filter(.isize >= min_overlap)
message(nrow(go_ds), " rows remain after requiring >= ", min_overlap,
        " overlapping genes")

## ---------------------------------------------------------------------------
## FOCAL / BARCHART TERMS
##   The same GO terms shown in GO_enrichment_OCD_CCD_facets.p5.barchart.pdf,
##   identified by applying the same filters as plot_fold_enrichment.R.
##   These are the ONLY focal terms (★) and the ONLY bolded rows.
## ---------------------------------------------------------------------------
barchart_terms <- go_ds %>%
  dplyr::filter(
    .comm  %in% communities,
    .pval  <  barchart_pvalue_max,
    .isize >  barchart_min_intersect
  ) %>%
  dplyr::group_by(.comm) %>%
  dplyr::arrange(dplyr::desc(.fe), .by_group = TRUE) %>%
  dplyr::slice_head(n = barchart_max_sets) %>%
  dplyr::ungroup() %>%
  dplyr::select(.comm, .go) %>%
  { split(.$.go, .$.comm) }

message("Focal (barchart) terms per community:")
for (cm in communities)
  message("  ", cm, ": ", length(barchart_terms[[cm]] %||% character(0)), " term(s)")

## ---------------------------------------------------------------------------
## GENE x COMMUNITY OVERLAP MATRIX
## ---------------------------------------------------------------------------
genes_by_comm <- setNames(
  lapply(communities, function(cm)
    sort(unique(unlist(go_ds$.genes[go_ds$.comm == cm])))),
  communities
)
all_overlap_genes <- sort(unique(unlist(genes_by_comm)))
if (length(all_overlap_genes) > 0) {
  overlap_mat <- tibble::tibble(gene = all_overlap_genes)
  for (cm in communities)
    overlap_mat[[cm]] <- overlap_mat$gene %in% genes_by_comm[[cm]]
  ord         <- order(-rowSums(as.matrix(overlap_mat[communities])), overlap_mat$gene)
  overlap_mat <- overlap_mat[ord, , drop = FALSE]
  readr::write_tsv(overlap_mat, paste0(plot_ds, "_community_gene_overlap.tsv"))
}

## ---------------------------------------------------------------------------
## SHARED HELPERS
## ---------------------------------------------------------------------------
get_safe <- function(id, map) tryCatch(AnnotationDbi::get(id, map),
                                       error = function(e) character(0))
san      <- function(x) str_replace_all(x, ":", "_")

# Abbreviate GO term names exactly as plot_fold_enrichment.v2.R does, plus
# "regulation" -> "reg." per user request, then wrap at the midpoint space.
abbreviate_go_label <- function(x) {
  x <- str_remove(x, " pathway$")                                           # strip suffix
  x <- str_replace_all(x, "G protein.coupled glutamate receptor", "mGluR") # must precede GPCR rule
  x <- str_replace_all(x, "G protein.coupled receptor",           "GPCR")  # saves 22 chars
  x <- str_replace_all(x, "adenylate cyclase",                    "AC")    # saves 15 chars
  x <- str_replace_all(x, "plasma membrane",                      "PM")    # saves 13 chars
  x <- str_replace_all(x, "[Cc]alcium.dependent",                 "Ca-dep.")
  x <- str_replace_all(x, "[Aa]dhesion",                          "adh.")
  x <- str_replace_all(x, "regulation of ",                       "reg. of ") # saves 6 chars
  x <- str_replace_all(x, "regulation",                           "reg.")  # remaining
  x <- paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))            # capitalize
  x
}

# No line-wrapping: just apply abbreviations.
wrap_name <- function(x) abbreviate_go_label(x)

classify_gene <- function(gene_u, net_syms) {
  is_net <- gene_u %in% net_syms
  is_dog <- gene_u %in% dog_seed_genes
  is_hum <- gene_u %in% human_seed_genes
  if      (is_net && is_dog) "network & dog seed"
  else if (is_net && is_hum) "network & human seed"
  else if (is_net)           "network gene"
  else if (is_dog)           "dog seed"
  else if (is_hum)           "human seed"
  else                       NA_character_
}

tree_layout <- function(tr) {
  ntip  <- length(tr$tip.label); nnode <- tr$Nnode; ntot <- ntip + nnode
  edge  <- tr$edge
  parent_of <- rep(NA_integer_, ntot); parent_of[edge[, 2]] <- edge[, 1]
  children  <- split(edge[, 2], edge[, 1])
  root <- unique(setdiff(edge[, 1], edge[, 2]))[1]
  depth <- rep(NA_real_, ntot); depth[root] <- 0; q <- root
  while (length(q)) {
    cur <- q[1]; q <- q[-1]
    ch  <- children[[as.character(cur)]]
    if (!is.null(ch)) { depth[ch] <- depth[cur] + 1; q <- c(q, ch) }
  }
  y <- rep(NA_real_, ntot); counter <- 0
  assign_y <- function(node) {
    ch <- children[[as.character(node)]]
    if (is.null(ch)) { counter <<- counter + 1; y[node] <<- counter; return(invisible()) }
    for (c in ch) assign_y(c)
    y[node] <<- mean(y[ch])
  }
  assign_y(root)
  label <- rep(NA_character_, ntot); label[seq_len(ntip)] <- tr$tip.label
  data.frame(node = seq_len(ntot), parent = parent_of, x = depth, y = y,
             isTip = seq_len(ntot) <= ntip, label = label,
             stringsAsFactors = FALSE)
}

## ---------------------------------------------------------------------------
## STEP 1: neighbourhood set for one focal term
## ---------------------------------------------------------------------------
compute_sel <- function(fg) {
  if (!(fg %in% sig_terms)) {
    warning("Focal term ", fg, " not in eligible pool - skipping.")
    return(NULL)
  }
  onto <- str_extract(sub$.source[sub$.go == fg][1], "BP|MF|CC")
  if (is.na(onto))
    onto <- tryCatch(as.character(AnnotationDbi::Ontology(GO.db::GOTERM[[fg]])),
                     error = function(e) NA_character_)
  if (is.na(onto) || !onto %in% c("BP", "MF", "CC")) {
    warning("Cannot resolve ontology for ", fg, " - skipping."); return(NULL)
  }
  ANCESTOR <- switch(onto, BP = GO.db::GOBPANCESTOR, MF = GO.db::GOMFANCESTOR, CC = GO.db::GOCCANCESTOR)
  PARENTS  <- switch(onto, BP = GO.db::GOBPPARENTS,  MF = GO.db::GOMFPARENTS,  CC = GO.db::GOCCPARENTS)
  CHILDREN <- switch(onto, BP = GO.db::GOBPCHILDREN, MF = GO.db::GOMFCHILDREN, CC = GO.db::GOCCCHILDREN)
  sel <- fg
  if (include_ancestors)
    sel <- union(sel, intersect(setdiff(get_safe(fg, ANCESTOR), "all"), sig_terms))
  else
    sel <- union(sel, intersect(get_safe(fg, PARENTS), sig_terms))
  if (include_siblings) {
    par  <- get_safe(fg, PARENTS)
    sibs <- unique(unlist(lapply(par, get_safe, map = CHILDREN)))
    sel  <- union(sel, intersect(sibs, sig_terms))
  }
  if (include_children)
    sel <- union(sel, intersect(get_safe(fg, CHILDREN), sig_terms))
  sel <- unique(c(fg, intersect(sel, sig_terms)))
  message("[", fg, "] ", onto, ": ", length(sel), " GO terms")
  list(focal = fg, onto = onto, sel = sel)
}

## ---------------------------------------------------------------------------
## STEP 2: build the GO-nesting tree
## ---------------------------------------------------------------------------
build_tree <- function(sel, onto) {
  ANCESTOR <- switch(onto, BP = GO.db::GOBPANCESTOR, MF = GO.db::GOMFANCESTOR, CC = GO.db::GOCCANCESTOR)
  in_set_anc <- setNames(
    lapply(sel, function(g) intersect(setdiff(get_safe(g, ANCESTOR), c("all", g)), sel)),
    sel)
  parent_of <- vapply(sel, function(g) {
    a <- in_set_anc[[g]]
    if (length(a) == 0) return(NA_character_)
    a[which.max(vapply(a, function(x) length(in_set_anc[[x]]), integer(1)))]
  }, character(1))
  names(parent_of) <- sel
  roots        <- sel[is.na(parent_of)]
  children_map <- setNames(lapply(sel, function(p)
    names(parent_of)[which(parent_of == p)]), sel)
  to_newick <- function(node) {
    kids <- children_map[[node]]
    if (length(kids) == 0) return(san(node))
    parts <- c(san(node), vapply(kids, to_newick, character(1)))
    paste0("(", paste(parts, collapse = ","), ")")
  }
  if (length(sel) == 1) {
    nwk <- paste0("(", san(sel), ");")
  } else if (length(roots) == 1) {
    nwk <- paste0(to_newick(roots), ";")
  } else {
    nwk <- paste0("(", paste(vapply(roots, to_newick, character(1)), collapse = ","), ");")
  }
  ape::read.tree(text = nwk)
}

## ---------------------------------------------------------------------------
## Build one panel.
##   focal_go (from enclosing scope): barchart terms — marked ★ AND bolded.
## ---------------------------------------------------------------------------
build_panel <- function(sel, tree, tile_x0, pad_to = NULL,
                        fe_limits = NULL, show_fe_label = TRUE,
                        panel_label = "", panel_height_in = NULL) {

  lay   <- tree_layout(tree)
  lay$x <- lay$x * tree_scale

  sel_sub    <- sub %>% dplyr::filter(.go %in% sel)
  net_by_go  <- setNames(sel_sub$.genes, sel_sub$.go)
  fe_by_go   <- setNames(sel_sub$.fe,    sel_sub$.go)
  name_by_go <- setNames(sel_sub$.name,  sel_sub$.go)
  pval_by_go <- setNames(sel_sub$.pval,  sel_sub$.go)
  netU_by_go <- lapply(net_by_go, norm_sym)

  seedextra_by_go <- setNames(lapply(sel, function(g) {
    ann   <- genes_in_term(g)
    setdiff(intersect(ann, c(dog_seed_genes, human_seed_genes)), netU_by_go[[g]])
  }), sel)

  disp_by_go <- setNames(lapply(sel, function(g)
    unique(c(net_by_go[[g]], seedextra_by_go[[g]]))), sel)

  all_genes  <- unique(unlist(disp_by_go))
  if (length(all_genes) == 0) all_genes <- character(0)
  gene_freq  <- table(factor(unlist(disp_by_go), levels = all_genes))
  net_all    <- unique(unlist(netU_by_go))
  is_net_col <- norm_sym(all_genes) %in% net_all
  gene_order <- all_genes[order(!is_net_col,
                                -as.integer(gene_freq[all_genes]),
                                all_genes)]
  n_net_cols <- sum(is_net_col)
  n_genes    <- length(gene_order)
  col_pos    <- seq_len(n_genes)
  if (n_net_cols > 0 && n_net_cols < n_genes)
    col_pos[(n_net_cols + 1):n_genes] <- col_pos[(n_net_cols + 1):n_genes] + 1
  col_x   <- setNames(col_pos, gene_order)
  n_slots <- if (n_genes) max(col_pos) else 0
  pad_to  <- pad_to %||% n_genes
  n_col   <- n_slots + max(0, pad_to - n_genes)

  tip_y <- setNames(lay$y[lay$isTip], lay$label[lay$isTip])

  tiledf <- do.call(rbind, lapply(sel, function(g) {
    gs <- intersect(disp_by_go[[g]], gene_order)
    if (!length(gs)) return(NULL)
    cats <- vapply(norm_sym(gs), classify_gene, character(1),
                   net_syms = netU_by_go[[g]])
    data.frame(x = tile_x0 + col_x[gs] - 0.5, y = tip_y[san(g)],
               fe = fe_by_go[[g]], go = g, gene = gs, category = cats,
               stringsAsFactors = FALSE)
  }))
  if (!is.null(tiledf)) {
    tiledf$category  <- factor(tiledf$category, levels = cat_levels)
    tiledf$dogcd_sig <- norm_sym(tiledf$gene) %in% dogcd_sig_genes
  }

  bg_grid <- expand.grid(go = sel, gene = gene_order,
                         KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  if (nrow(bg_grid)) {
    bg_grid$x <- tile_x0 + col_x[bg_grid$gene] - 0.5
    bg_grid$y <- unname(tip_y[san(bg_grid$go)])
    bg_grid   <- bg_grid[!is.na(bg_grid$y), , drop = FALSE]
  }

  ed   <- lay[!is.na(lay$parent), ]
  px   <- lay$x[match(ed$parent, lay$node)]
  hseg <- data.frame(x = px, xend = ed$x, y = ed$y, yend = ed$y)
  vseg <- do.call(rbind, lapply(unique(ed$parent), function(pp) {
    ys <- ed$y[ed$parent == pp]
    data.frame(x = lay$x[lay$node == pp], xend = lay$x[lay$node == pp],
               y = min(ys), yend = max(ys))
  }))

  # ★ marker and bold face: only for focal (barchart) terms
  star    <- ifelse(sel %in% focal_go, "★ ", "")
  disp    <- setNames(
    paste0(star, sel, "  ", wrap_name(name_by_go[sel]),
           "  (p=", formatC(pval_by_go[sel], format = "e", digits = 1), ")"),
    san(sel))
  tip_lab <- lay[lay$isTip, ]
  tip_lab <- tip_lab[order(tip_lab$y), ]
  y_face  <- ifelse(tip_lab$label %in% san(focal_go), "bold", "plain")

  x_breaks <- tile_x0 + seq_len(n_col) - 0.5
  x_labels <- rep("", n_col); x_labels[col_pos] <- gene_order

  seed_cols <- if (n_net_cols < n_genes) gene_order[(n_net_cols + 1):n_genes] else character(0)
  box_df <- NULL
  if (length(seed_cols) > 0 && length(tip_y) > 0) {
    s_slots <- col_x[seed_cols]
    box_df <- data.frame(
      xmin = tile_x0 + min(s_slots) - 1 - 0.18,
      xmax = tile_x0 + max(s_slots) + 0.18,
      ymin = min(tip_y) - 0.65, ymax = max(tip_y) + 0.65)
  }

  fe_axis_max <- if (!is.null(fe_limits) && is.finite(max(fe_limits, na.rm = TRUE)))
    max(fe_limits, na.rm = TRUE) else suppressWarnings(max(fe_by_go[sel], na.rm = TRUE))
  if (!is.finite(fe_axis_max) || fe_axis_max <= 0) fe_axis_max <- 1
  bar_x0 <- tile_x0 - bar_gap - bar_w_units

  # Bar half-height in data coordinates, normalised so bars have the same
  # physical height across all panels.
  #
  # The displayed y range is n_tips + 2.55 data units (derived from:
  #   expand = expansion(add = c(1.5, 0.5)), tick labels at min(tip_y)-1.05,
  #   bar-axis segment top at max(tip_y)+0.5).
  # We target bar_half = 0.45 (same as geom_tile height/2 = 0.9/2).
  # Physical bar height = 2 * bar_half * (data_area_in / data_range)
  #   where data_area_in ≈ panel_height_in - axis_text_oh.
  # Substituting the formula for bar_half, all terms cancel to give
  # physical height = 0.9 * row_height_in * height_scale (constant).
  n_tips   <- length(sel)
  # With heights_in = (n_rows + y_data_pad) * row_height_in * height_scale + axis_oh_in,
  # the data area = (n_tips + y_data_pad) * row_height_in * height_scale, so
  # bar_half = 0.45 * data_range / data_area cancels to exactly 0.45 (= tile height/2).
  bar_half <- if (!is.null(panel_height_in) && panel_height_in > 0) {
    eff_h <- max(panel_height_in - axis_oh_in,
                 n_tips * row_height_in * height_scale * 0.5)
    0.45 * row_height_in * height_scale * (n_tips + y_data_pad) / eff_h
  } else
    0.45

  bardf <- data.frame(y = unname(tip_y[san(sel)]), fe = unname(fe_by_go[sel]),
                      is_focal = sel %in% focal_go,
                      stringsAsFactors = FALSE)
  bardf <- bardf[!is.na(bardf$y) & !is.na(bardf$fe), , drop = FALSE]
  bardf$xmin <- bar_x0
  bardf$xmax <- bar_x0 + (bardf$fe / fe_axis_max) * bar_w_units
  # Bake bar_half into the data frame so ggplot never evaluates it lazily from
  # an outer environment (which can resolve to the wrong value at render time).
  bardf$ymin <- bardf$y - bar_half
  bardf$ymax <- bardf$y + bar_half

  fe_breaks <- pretty(c(0, fe_axis_max), n = 4)
  fe_breaks <- fe_breaks[fe_breaks >= 0 & fe_breaks <= fe_axis_max + 1e-9]
  bar_y_top <- if (length(tip_y)) max(tip_y) + 0.5 else 1
  bar_y_bot <- if (length(tip_y)) min(tip_y) - 0.5 else 0
  axis_y    <- bar_y_bot - 0.55
  title_y   <- bar_y_bot - 1.15
  gx        <- bar_x0 + (fe_breaks / fe_axis_max) * bar_w_units
  grid_df   <- data.frame(x = gx, y = bar_y_bot, yend = bar_y_top)
  tick_df   <- data.frame(x = gx, y = axis_y, lab = formatC(fe_breaks, format = "g"))

  long <- if (is.null(tiledf)) {
    tibble::tibble(panel = character(0), go = character(0), gene = character(0),
                   fe = numeric(0), category = character(0))
  } else {
    tibble::tibble(panel = panel_label, go = tiledf$go, gene = tiledf$gene,
                   fe = tiledf$fe, category = as.character(tiledf$category))
  }

  gg <- ggplot() +
    geom_segment(data = hseg, aes(x = x, y = y, xend = xend, yend = yend),
                 colour = "grey45", linewidth = 0.3) +
    geom_segment(data = vseg, aes(x = x, y = y, xend = xend, yend = yend),
                 colour = "grey45", linewidth = 0.3) +
    geom_tile(data = bg_grid, aes(x = x, y = y), width = 1, height = 0.9,
              fill = NA, colour = "grey90", linewidth = 0.2, inherit.aes = FALSE) +
    geom_tile(data = tiledf, aes(x = x, y = y, fill = category),
              width = 1, height = 0.9, colour = "grey92", linewidth = 0.2) +
    { if (!is.null(tiledf) && any(tiledf$dogcd_sig))
        geom_text(data = tiledf[tiledf$dogcd_sig, , drop = FALSE],
                  aes(x = x, y = y), label = "*",
                  size = 2.8, colour = "white", fontface = "bold",
                  vjust = 0.7, alpha = 0.8, inherit.aes = FALSE)
      else NULL } +
    geom_segment(data = grid_df, aes(x = x, y = y, xend = x, yend = yend),
                 colour = "grey88", linewidth = 0.3, inherit.aes = FALSE) +
    geom_rect(data = bardf[!bardf$is_focal, , drop = FALSE],
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = "grey65", colour = NA, inherit.aes = FALSE) +
    geom_rect(data = bardf[bardf$is_focal, , drop = FALSE],
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = "black", colour = NA, inherit.aes = FALSE) +
    geom_text(data = tick_df, aes(x = x, y = y, label = lab),
              size = 1.9, colour = "grey40", vjust = 1, inherit.aes = FALSE) +
    { if (show_fe_label) annotate("text", x = bar_x0 + bar_w_units / 2, y = title_y,
             label = "Fold enrichment", size = 2.4, colour = "grey30", vjust = 1) else NULL } +
    scale_fill_manual(name = "Gene type", values = cat_colors,
                      breaks = cat_levels, drop = FALSE, na.value = "white") +
    scale_x_continuous(position = "top", breaks = x_breaks, labels = x_labels,
                       expand = expansion(mult = 0)) +
    scale_y_continuous(breaks = tip_lab$y, labels = disp[tip_lab$label],
                       expand = expansion(add = c(1.5, 0.5))) +
    coord_cartesian(xlim = c(-0.3, tile_x0 + n_col + 0.3), clip = "off") +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 9) +
    theme(
      axis.text.x  = element_text(angle = 90, hjust = 0, vjust = 0.5, size = 6),
      axis.text.y  = element_text(size = 7, face = y_face),
      panel.grid   = element_blank(),
      legend.position = "right",
      plot.margin  = margin(0, 0, 0, 0)
    )

  if (!is.null(box_df))
    gg <- gg + geom_rect(data = box_df,
                         aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                         inherit.aes = FALSE, fill = NA,
                         colour = "grey62", linewidth = 0.6)

  list(gg = gg, n = length(sel), long = long)
}

## ---------------------------------------------------------------------------
## PASS 1 — collect per-community layout data
## ---------------------------------------------------------------------------
comm_data <- list()

for (cm in communities) {

  sub <- go_ds %>%
    dplyr::filter(.comm == cm) %>%
    dplyr::group_by(.go) %>% dplyr::slice(1) %>% dplyr::ungroup()
  if (nrow(sub) == 0) { message("== ", cm, ": no eligible rows - skipping"); next }

  sig_terms <- sub$.go
  fe_of     <- setNames(sub$.fe, sub$.go)

  # focal terms = barchart terms for this community (intersected with eligible pool)
  focal_go <- intersect(barchart_terms[[cm]] %||% character(0), sig_terms)
  if (length(focal_go) == 0) {
    message("== ", cm, ": no barchart focal terms in eligible pool - skipping"); next
  }
  message("== ", cm, ": ", length(focal_go), " focal term(s) ==")

  sels <- Filter(Negate(is.null), lapply(focal_go, compute_sel))
  if (length(sels) == 0) { message("  no resolvable focals for ", cm); next }
  nfoc <- length(sels)

  if (merge_focals) {
    sel_list <- lapply(sels, `[[`, "sel")
    onto_v   <- vapply(sels, `[[`, character(1), "onto")
    foc_v    <- vapply(sels, `[[`, character(1), "focal")

    anc_depths <- function(go, onto, maxstep) {
      PARENTS <- switch(onto, BP = GO.db::GOBPPARENTS, MF = GO.db::GOMFPARENTS, CC = GO.db::GOCCPARENTS)
      depth <- setNames(0L, go); frontier <- go; step <- 0L
      while (step < maxstep && length(frontier)) {
        step <- step + 1L
        nxt <- setdiff(unique(unlist(lapply(frontier, get_safe, map = PARENTS))),
                       c("all", names(depth)))
        if (!length(nxt)) break
        depth <- c(depth, setNames(rep(step, length(nxt)), nxt))
        frontier <- nxt
      }
      depth
    }
    depl <- lapply(seq_len(nfoc), function(i) anc_depths(foc_v[i], onto_v[i], merge_dist))

    cand <- list()
    if (nfoc > 1)
      for (i in 1:(nfoc - 1)) for (j in (i + 1):nfoc) {
        if (onto_v[i] != onto_v[j]) next
        common <- intersect(names(depl[[i]]), names(depl[[j]]))
        if (!length(common)) next
        d <- min(depl[[i]][common] + depl[[j]][common])
        if (is.finite(d) && d <= merge_dist) cand[[length(cand) + 1]] <- c(i, j, d)
      }

    parent  <- seq_len(nfoc)
    froot   <- function(i) { while (parent[i] != i) { parent[i] <<- parent[parent[i]]; i <- parent[i] }; i }
    grp_sel <- sel_list
    if (length(cand)) {
      cand <- cand[order(vapply(cand, `[[`, numeric(1), 3))]
      for (lk in cand) {
        ri <- froot(lk[1]); rj <- froot(lk[2])
        if (ri == rj) next
        merged <- union(grp_sel[[ri]], grp_sel[[rj]])
        if (length(merged) <= max_panel_terms) {
          parent[ri] <- rj; grp_sel[[rj]] <- merged
        }
      }
    }
    groups <- split(seq_len(nfoc), vapply(seq_len(nfoc), froot, integer(1)))
    keep <- lapply(unname(groups), function(ix) {
      foci <- foc_v[ix]
      sel  <- intersect(unique(unlist(sel_list[ix])), sig_terms)
      list(sel   = sel, onto  = onto_v[ix[1]], foci  = foci,
           genes = unique(unlist(sub$.genes[sub$.go %in% sel])),
           maxfe = suppressWarnings(max(fe_of[foci], na.rm = TRUE)),
           label = foci[which.max(fe_of[foci])])
    })
    keep <- keep[order(-vapply(keep, `[[`, numeric(1), "maxfe"))]
    message("  ", nfoc, " focal(s) -> ", length(keep), " panel(s)")
  } else {
    sels    <- sels[order(-fe_of[vapply(sels, `[[`, character(1), "focal")])]
    covered <- character(0); keep <- list()
    for (s in sels) {
      if (length(intersect(s$sel, covered)) > 0) next
      s$foci <- s$focal; s$label <- s$focal
      s$genes <- unique(unlist(sub$.genes[sub$.go %in% s$sel]))
      keep[[length(keep) + 1]] <- s; covered <- union(covered, s$sel)
    }
  }

  panel_width <- function(s) {
    net  <- unique(unlist(sub$.genes[sub$.go %in% s$sel]))
    netU <- norm_sym(net)
    extra <- setdiff(intersect(unique(unlist(lapply(s$sel, genes_in_term))),
                               c(dog_seed_genes, human_seed_genes)), netU)
    length(unique(c(net, extra)))
  }

  trees        <- lapply(keep, function(s) build_tree(s$sel, s$onto))
  depths       <- vapply(trees, function(t) max(tree_layout(t)$x), numeric(1))
  tile_x0_cm   <- max(depths) * tree_scale + tree_gap + bar_w_units + bar_gap
  max_genes_cm <- max(vapply(keep, panel_width, numeric(1)))
  fe_limits_cm <- range(fe_of[unlist(lapply(keep, `[[`, "sel"))], na.rm = TRUE)

  comm_data[[cm]] <- list(
    sub       = sub, sig_terms = sig_terms, focal_go  = focal_go,
    fe_of     = fe_of, keep = keep, trees = trees,
    tile_x0   = tile_x0_cm, max_genes = max_genes_cm, fe_limits = fe_limits_cm
  )
}

if (length(comm_data) == 0) stop("No communities produced any panels.")

## ---------------------------------------------------------------------------
## GLOBAL LAYOUT PARAMETERS
## ---------------------------------------------------------------------------
global_tile_x0   <- max(vapply(comm_data, `[[`, numeric(1), "tile_x0"))
global_max_genes <- max(vapply(comm_data, `[[`, numeric(1), "max_genes"))
global_fe_limits <- range(unlist(lapply(comm_data, function(cd) cd$fe_limits)),
                          na.rm = TRUE)
message("Global tile_x0=", round(global_tile_x0, 2),
        "  max_genes=", global_max_genes,
        "  fe_range=", paste(round(global_fe_limits, 1), collapse = "-"))

## ---------------------------------------------------------------------------
## PASS 2 — build panels with globally aligned layout
## ---------------------------------------------------------------------------
comm_ggs     <- list()
comm_heights <- list()
all_long     <- list()

for (cm in names(comm_data)) {
  cd          <- comm_data[[cm]]
  sub         <- cd$sub
  sig_terms   <- cd$sig_terms
  focal_go    <- cd$focal_go
  fe_of       <- cd$fe_of
  n_panels_cm <- length(cd$keep)

  # Pre-compute heights so panel_height_in can be passed to build_panel.
  # Formula: (n_rows + y_data_pad) * row_height_in * height_scale + axis_oh_in
  # This makes every panel have the same physical row height (inches per data unit
  # is constant), so tiles and bars are visually consistent across all panels.
  n_rows_v   <- vapply(cd$keep, function(s) length(s$sel), integer(1))
  heights_in <- (n_rows_v + y_data_pad) * row_height_in * height_scale + axis_oh_in

  panels <- Map(function(s, tr, idx, h) {
    build_panel(s$sel, tree = tr,
                tile_x0         = global_tile_x0,
                pad_to          = global_max_genes,
                fe_limits       = global_fe_limits,
                show_fe_label   = (idx == n_panels_cm),
                panel_label     = paste0("Community ", cm, "  #", idx),
                panel_height_in = h)
  }, cd$keep, cd$trees, seq_along(cd$keep), heights_in)

  n_rows <- n_rows_v

  comm_ggs[[cm]]     <- lapply(panels, `[[`, "gg")
  comm_heights[[cm]] <- heights_in
  all_long           <- c(all_long, lapply(panels, `[[`, "long"))
  message("  ", cm, ": ", length(panels), " panel(s), ", sum(n_rows), " rows")
}

## ---------------------------------------------------------------------------
## COMBINED ONE-PAGE OUTPUT
## ---------------------------------------------------------------------------
legend_src <- ggplot(
  data.frame(cat = factor(cat_levels, levels = cat_levels),
             x = 1, y = seq_along(cat_levels)),
  aes(x, y, fill = cat)) +
  geom_tile() +
  scale_fill_manual(name = "Gene type", values = cat_colors,
                    breaks = cat_levels, drop = FALSE,
                    guide = guide_legend(nrow = 1)) +
  theme_minimal(base_size = 9) +
  theme(legend.position   = "top",  legend.direction  = "horizontal",
        legend.box.margin = margin(0, 0, 0, 0),
        legend.key.size   = grid::unit(0.85, "lines"),
        legend.title      = element_text(size = 8),
        legend.text       = element_text(size = 7))
legend_grob <- cowplot::get_legend(legend_src)

strip_legend <- function(g) g + theme(legend.position = "none")

# Flatten all panels in community order so align="v" acts across ALL communities
all_panels_flat  <- unlist(lapply(names(comm_ggs),
                                  function(cm) lapply(comm_ggs[[cm]], strip_legend)),
                           recursive = FALSE)
all_heights_flat <- unlist(comm_heights, use.names = FALSE)

# Body column: single plot_grid so horizontal spacing is globally aligned
body_col <- cowplot::plot_grid(
  plotlist    = all_panels_flat,
  ncol        = 1,
  rel_heights = all_heights_flat,
  align       = "v",
  axis        = "lr"
)

# Label column: one rotated label per community in a coloured box
comm_total_heights <- vapply(names(comm_ggs),
                             function(cm) sum(comm_heights[[cm]]), numeric(1))
make_comm_label <- function(cm) {
  fill <- comm_colors[[cm]] %||% "grey40"
  cowplot::ggdraw() +
    cowplot::draw_grob(grid::rectGrob(gp = grid::gpar(fill = fill, col = NA))) +
    cowplot::draw_label(cm, x = 0.5, y = 0.5, angle = 90,
                        fontface = "bold", size = 10, color = "white")
}
label_col <- cowplot::plot_grid(
  plotlist    = lapply(names(comm_ggs), make_comm_label),
  ncol        = 1,
  rel_heights = comm_total_heights
)

body <- cowplot::plot_grid(label_col, body_col, nrow = 1,
                           rel_widths = c(0.025, 0.975))

all_heights_in <- comm_total_heights
legend_h <- 0.5
combined <- cowplot::plot_grid(legend_grob, body, ncol = 1,
                               rel_heights = c(legend_h, sum(all_heights_in)))
fig_h    <- sum(all_heights_in) + legend_h

all_terms  <- unique(unlist(lapply(comm_data, function(cd)
  unique(unlist(lapply(cd$keep, `[[`, "sel"))))))
all_names  <- go_ds$.name[match(all_terms, go_ds$.go)]
all_names[is.na(all_names)] <- ""
ylab_chars <- max(nchar(all_terms) + nchar(all_names) + 19)
ylab_in    <- ylab_chars * ylab_char_in + 0.15
plot_w_use <- if (square_tiles) {
  (global_tile_x0 + global_max_genes * 0.9 + 1.6) * row_height_in * height_scale + ylab_in
} else {
  14
}

out_stem <- paste0(plot_ds, "_all_communities_p4_maxts", max_termsize)
out_pdf  <- paste0(out_stem, "_nesting_tileplot.pdf")
out_csv  <- paste0(out_stem, "_nesting_tileplot_data.csv")

ggsave(out_pdf, plot = combined, width = plot_w_use,
       height = fig_h, limitsize = FALSE)
readr::write_csv(
  dplyr::bind_rows(all_long) %>%
    dplyr::transmute(community = str_extract(panel, "C\\d+"),
                     panel, go, gene, fold_enrichment = fe, category),
  out_csv
)

message("Done. Wrote ", out_pdf, " (", sum(lengths(comm_ggs)), " panel(s) on one page)")
message("      Data: ", out_csv)

## ---------------------------------------------------------------------------
## LEGEND INFO: abbreviations actually used
## ---------------------------------------------------------------------------
abbrev_glossary <- c(
  "mGluR"   = "G protein-coupled glutamate receptor",
  "GPCR"    = "G protein-coupled receptor",
  "AC"      = "adenylate cyclase",
  "PM"      = "plasma membrane",
  "Ca-dep." = "calcium-dependent",
  "adh."    = "adhesion",
  "reg."    = "regulation"
)

shown_names  <- unique(go_ds$.name[go_ds$.go %in% all_terms & !is.na(go_ds$.name)])
orig_abbr    <- vapply(shown_names, abbreviate_go_label, character(1))
changed_mask <- shown_names != orig_abbr
used_pairs   <- data.frame(
  original    = shown_names[changed_mask],
  abbreviated = orig_abbr[changed_mask],
  stringsAsFactors = FALSE
)
# Determine which glossary entries were responsible
used_pairs$rule <- vapply(seq_len(nrow(used_pairs)), function(i) {
  hits <- names(abbrev_glossary)[
    sapply(abbrev_glossary, function(full) grepl(full, used_pairs$original[i], ignore.case = TRUE))
  ]
  if (length(hits)) paste(hits, collapse = ", ") else "other"
}, character(1))

out_legend <- sub("\\.pdf$", "_legend.info.txt", out_pdf)
writeLines(c(
  paste0("Abbreviation legend for: ", out_pdf),
  paste0("Generated: ", Sys.time()),
  "",
  "--- Glossary (all defined abbreviations) ---",
  paste0(formatC(names(abbrev_glossary), width = 8, flag = "-"),
         " = ", unname(abbrev_glossary)),
  "",
  "--- Terms abbreviated in this figure ---",
  if (nrow(used_pairs) > 0)
    paste0("  ", used_pairs$original, "  ->  ", used_pairs$abbreviated)
  else
    "  (none)"
), out_legend)
message("      Legend: ", out_legend)
