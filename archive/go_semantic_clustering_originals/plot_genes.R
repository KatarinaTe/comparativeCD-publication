#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# go_nesting_tileplot.R
#
# Tile plot of GO sets (rows) x genes (columns) with a GO-hierarchy tree on the
# y-axis, focused on a target GO term plus its significant ancestors/siblings.
# Tile fill = fold enrichment of the GO set the gene belongs to.
#
# Reads gProfiler-style combined results from cache/all_go_combined.rds.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(ggtree)   # Bioconductor
  library(GO.db)    # Bioconductor
  library(ape)
  library(aplot)
})
# NB: GO.db loads AnnotationDbi, which masks dplyr::select / dplyr::filter.
#     All dplyr verbs below are namespaced (dplyr::) to avoid that collision.

`%||%` <- function(a, b) if (is.null(a)) b else a

## ---------------------------------------------------------------------------
## CONFIG  -- edit these
## ---------------------------------------------------------------------------
cache_path <- "all_go_combined.rds"
plot_ds    <- "OCD_CCD"
community  <- "C185"
focal_go   <- "GO:1904862"   # ionotropic glutamate receptor signalling pathway

p_cutoff   <- 0.01           # only include GO sets with p_value < this
max_termsize <- 2000         # exclude GO sets larger than this many genes

include_ancestors <- TRUE    # full significant lineage above focal (FALSE = direct parents only)
include_siblings  <- TRUE    # significant terms sharing a parent with focal
include_children  <- FALSE   # significant direct children of focal

out_stem   <- paste0(plot_ds, "_", community, "_", gsub(":", "", focal_go),
                     "_p", -log10(p_cutoff), "_max", max_termsize)
out_pdf    <- paste0(out_stem, "_nesting_tileplot.pdf")
out_csv    <- paste0(out_stem, "_nesting_tileplot_data.csv")
plot_w     <- 16
plot_h     <- 8
tree_width <- 0.18           # tree width as a fraction of the combined plot

# Column names. NULL entries are auto-detected; glimpse() is printed at the
# start so you can hard-set any that the detector gets wrong.
col_ds       <- NULL         # NULL -> column whose values include plot_ds
col_comm     <- "community"
col_go       <- "native"
col_name     <- "name"
col_sig      <- "significant"
col_pval     <- "p_value"
col_termsize <- "term_size"
col_qsize    <- "query_size"
col_source   <- "source"
col_genes    <- NULL         # NULL -> the bracketed gene-list column
col_effdom   <- NULL         # NULL -> the per-row effective_domain_size column

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
  if (length(hit) > 1) {                      # prefer gene-symbol-looking tokens
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

sub <- go %>%
  dplyr::filter(.data[[col_ds]]       == plot_ds,
                .data[[col_comm]]     == community,
                .data[[col_sig]]      == TRUE,
                .data[[col_pval]]     <  p_cutoff,
                .data[[col_termsize]] <= max_termsize)
message(nrow(sub), " significant rows (p < ", p_cutoff,
        ", term_size <= ", max_termsize, ") for ds=",
        plot_ds, " community=", community)
stopifnot(nrow(sub) > 0)

col_effdom <- col_effdom %||%
  detect_effdom_col(sub, exclude = c(col_termsize, col_qsize))
if (is.na(col_effdom))
  stop("Could not find an effective_domain_size column; set col_effdom manually.")
message("effective_domain_size column: ", col_effdom)

sub <- sub %>%
  dplyr::mutate(
    .go     = .data[[col_go]],
    .name   = .data[[col_name]],
    .source = .data[[col_source]],
    .genes  = parse_genes(.data[[col_genes]]),
    .isize  = lengths(.genes),
    # fold enrichment = (k/n) / (K/N)
    .fe     = (.isize / .data[[col_qsize]]) /
      (.data[[col_termsize]] / .data[[col_effdom]])
  ) %>%
  dplyr::group_by(.go) %>% dplyr::slice(1) %>% dplyr::ungroup()

stopifnot(focal_go %in% sub$.go)   # focal must be significant + present here

## ---------------------------------------------------------------------------
## GO NEIGHBOURHOOD (focal + significant ancestors / siblings / children)
## ---------------------------------------------------------------------------
onto <- str_extract(sub$.source[sub$.go == focal_go][1], "BP|MF|CC")
if (is.na(onto)) onto <- as.character(Ontology(GOTERM[[focal_go]]))
message("Focal ", focal_go, " ontology: ", onto)

ANCESTOR <- switch(onto, BP = GOBPANCESTOR, MF = GOMFANCESTOR, CC = GOCCANCESTOR)
PARENTS  <- switch(onto, BP = GOBPPARENTS,  MF = GOMFPARENTS,  CC = GOCCPARENTS)
CHILDREN <- switch(onto, BP = GOBPCHILDREN, MF = GOMFCHILDREN, CC = GOCCCHILDREN)

sig_terms <- sub$.go
get_safe  <- function(id, map) tryCatch(AnnotationDbi::get(id, map),
                                        error = function(e) character(0))

sel <- focal_go
if (include_ancestors) {
  sel <- union(sel, intersect(setdiff(get_safe(focal_go, ANCESTOR), "all"), sig_terms))
} else {
  sel <- union(sel, intersect(get_safe(focal_go, PARENTS), sig_terms))
}
if (include_siblings) {
  par  <- get_safe(focal_go, PARENTS)
  sibs <- unique(unlist(lapply(par, get_safe, map = CHILDREN)))
  sel  <- union(sel, intersect(sibs, sig_terms))
}
if (include_children) {
  sel <- union(sel, intersect(get_safe(focal_go, CHILDREN), sig_terms))
}
sel <- unique(c(focal_go, intersect(sel, sig_terms)))
message("Selected ", length(sel), " GO terms:\n  ", paste(sel, collapse = "\n  "))

## ---------------------------------------------------------------------------
## REDUCE THE GO DAG TO A TREE
##   - each term is attached to its NEAREST in-set ancestor
##   - every term (including internal ancestors) becomes a tip via a self-leaf,
##     so each GO set gets its own row in the heatmap
##   - ":" is illegal in Newick, so labels are sanitised and mapped back
## ---------------------------------------------------------------------------
in_set_anc <- setNames(
  lapply(sel, function(g) intersect(setdiff(get_safe(g, ANCESTOR), c("all", g)), sel)),
  sel
)
parent_of <- vapply(sel, function(g) {
  a <- in_set_anc[[g]]
  if (length(a) == 0) return(NA_character_)
  a[which.max(vapply(a, function(x) length(in_set_anc[[x]]), integer(1)))]
}, character(1))
names(parent_of) <- sel

roots        <- sel[is.na(parent_of)]
children_map <- setNames(lapply(sel, function(p)
  names(parent_of)[which(parent_of == p)]), sel)

san       <- function(x) str_replace_all(x, ":", "_")
to_newick <- function(node) {
  kids <- children_map[[node]]
  if (length(kids) == 0) return(san(node))                 # pure leaf
  parts <- c(san(node), vapply(kids, to_newick, character(1)))  # self-leaf + kids
  paste0("(", paste(parts, collapse = ","), ")")           # unlabelled internal node
}
nwk  <- if (length(roots) == 1) {
  paste0(to_newick(roots), ";")
} else {
  paste0("(", paste(vapply(roots, to_newick, character(1)), collapse = ","), ");")
}
tree <- ape::read.tree(text = nwk)

## ---------------------------------------------------------------------------
## ASSEMBLE PLOT DATA
## ---------------------------------------------------------------------------
sel_sub    <- sub %>% dplyr::filter(.go %in% sel)
sel_genes  <- setNames(sel_sub$.genes, sel_sub$.go)
fe_by_go   <- setNames(sel_sub$.fe,    sel_sub$.go)
name_by_go <- setNames(sel_sub$.name,  sel_sub$.go)

all_genes  <- unique(unlist(sel_genes))
gene_freq  <- table(factor(unlist(sel_genes), levels = all_genes))
gene_order <- all_genes[order(-as.integer(gene_freq[all_genes]), all_genes)]

long <- imap_dfr(sel_genes, ~ tibble(go = .y, gene = .x)) %>%
  dplyr::mutate(
    fe   = fe_by_go[go],
    san  = san(go),
    gene = factor(gene, levels = gene_order)
  )

# sanitised-id -> display label ("GO:xxxxxxx  term name") for the y axis
disp_lab <- setNames(paste0(sel, "  ", name_by_go[sel]), san(sel))

## ---------------------------------------------------------------------------
## PLOT  (tree | GO labels | tiles)
## ---------------------------------------------------------------------------
tile <- ggplot(long, aes(x = gene, y = san, fill = fe)) +
  geom_tile(colour = "grey92", linewidth = 0.2) +
  scale_fill_viridis_c(name = "Fold\nenrichment", option = "C", na.value = "white") +
  scale_y_discrete(labels = disp_lab) +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL,
       title = paste0(plot_ds, " ", community, " - ", focal_go,
                      " neighbourhood")) +
  theme_minimal(base_size = 9) +
  theme(
    axis.text.x      = element_text(angle = 90, hjust = 0, vjust = 0.5, size = 6),
    axis.text.y      = element_text(size = 7),
    panel.grid       = element_blank(),
    legend.position  = "right",
    plot.title       = element_text(size = 10, face = "bold")
  )

tree_plot <- ggtree(tree)

final <- tile %>% insert_left(tree_plot, width = tree_width)

ggsave(out_pdf, plot = final, width = plot_w, height = plot_h, limitsize = FALSE)
readr::write_csv(
  long %>% dplyr::transmute(go, term = name_by_go[go], gene, fold_enrichment = fe),
  out_csv
)
message("Wrote ", out_pdf, " and ", out_csv)