#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# go_nesting_tileplot.R
#
# Tile plot of GO sets (rows) x genes (columns) with a GO-hierarchy tree on the
# y-axis, focused on a target GO term plus its significant ancestors/siblings.
#
# Tile colour (fill) = gene category:
#     "network gene", "dog seed", "human seed",
#     "network & dog seed", "network & human seed"
# A horizontal bar to the left of each row shows that GO set's fold enrichment,
# on a common axis shared across all panels in a PDF.
#
# Seed genes (dog / human) are read from a frozen CSV export of a Google Sheet. A seed gene is shown
# on the plot if it is directly annotated to a plotted GO term in the GO
# Consortium human GAF, even if it was not one of the network (query-
# intersection) genes.
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
  library(ape)
  library(cowplot)
})
# GO.db / AnnotationDbi supply the GO hierarchy (ancestor / parent / child maps
# and a term's ontology). They are referenced with the `pkg::object` form (e.g.
# GO.db::GOBPANCESTOR, AnnotationDbi::Ontology), which loads the namespace on
# first use without attaching it to the search path, so dplyr::select /
# dplyr::filter stay unmasked. Both packages must be installed.
# Gene -> GO annotations come from the downloadable GO Consortium human GAF
# (see goa_gaf below); org.Hs.eg.db is not used.
stopifnot(requireNamespace("GO.db", quietly = TRUE),
          requireNamespace("AnnotationDbi", quietly = TRUE))

`%||%` <- function(a, b) if (is.null(a)) b else a

## ---------------------------------------------------------------------------
## CONFIG  -- edit these
## ---------------------------------------------------------------------------
cache_path <- "all_go_combined.rds"
plot_ds    <- "OCD_CCD"

# --- seed genes (frozen CSV) ------------------------------------------------
# The "gene list w gwas catalog overlap" tab provides per-gene flags. Seed sets
# are taken from the logical columns `dog seed` and `human seed`; gene symbols
# come from the `gene` column.
seed_file   <- "gene_list_gwas_catalog_overlap_with_region.csv"  # frozen export of the seed-gene sheet tab below
seed_tab    <- "gene list w gwas catalog overlap"
col_gene_id <- "gene"        # column holding the gene symbol
col_dogseed <- "dog seed"    # logical column: TRUE = dog seed gene
col_humseed <- "human seed"  # logical column: TRUE = human seed gene
# Dog seed symbols are assumed to already be human-style symbols and are matched
# directly (case-insensitively) against network / GO gene symbols.

# --- authoritative GO annotations (GO Consortium human GAF) ------------------
# Seed-gene -> GO-term membership uses the GO Consortium's human GAF: a gene
# belongs to a term only if it is directly annotated to that exact term (no
# propagation to descendant terms). The file is downloaded once and cached
# locally; delete it (or bump the URL) to refresh. To use a file you've already
# downloaded, just set goa_gaf to its path.
goa_gaf <- "goa_human.gaf.gz"
goa_url <- "http://current.geneontology.org/annotations/goa_human.gaf.gz"

# One PDF is produced per community. Focal GO terms are picked automatically
# per community (see "interesting term" thresholds below) rather than supplied.
# Set to a character vector to choose communities explicitly, e.g.
#   communities <- c("C185", "C186", "C197", "C187", "C193", "C403")
# Leave as NULL (or empty) to auto-pick the 5 biggest clusters (see below).
communities <- c("C184","C185", "C186", "C193", "C197")
n_biggest    <- max(length(communities),5)            # how many clusters to use when `communities` is NULL

# --- "interesting" focal-term selection (per community) ---------------------
pvalue_max       <- 1e-5     # focal terms must have p_value < this
min_intersection <- 5        # ... and intersection_size > this
min_term_size    <- 5        # ... and term_size > this
# (no cap on the number of focal terms: every term passing the criteria is used)

# --- pool of terms eligible to appear as rows (focals + their relatives) ----
p_cutoff     <- 1e-4         # ancestors/siblings shown must have p_value < this
max_termsize <- 500         # exclude GO sets larger than this many genes
min_overlap  <- 3            # exclude any term with fewer than this many overlapping genes

include_ancestors <- TRUE    # full significant lineage above focal (FALSE = direct parents only)
include_siblings  <- TRUE    # significant terms sharing a parent with focal
include_children  <- FALSE   # significant direct children of focal

# --- merging focal terms into fewer panels per PDF --------------------------
# When merge_focals is TRUE, focal terms that sit close together in the GO
# hierarchy are grouped into a shared panel (one tree), so related focals are
# not drawn as separate panels and none are dropped. Two dials control it:
#   merge_dist      -- how close two focals must be to be linked, measured as
#                      the number of GO steps from each focal up to their nearest
#                      shared ancestor, summed over both sides (parent<->child=1,
#                      two siblings under one parent=2, ...). Larger = more
#                      grouping = fewer panels.
#   max_panel_terms -- hard cap on the number of GO-term rows in one panel. A
#                      cluster larger than this is split into several panels
#                      (closest focals merged first), so no single panel becomes
#                      huge. Lower = smaller, more numerous panels.
# Set merge_focals to FALSE for the original behaviour (one panel per focal,
# greedily dropping any focal already covered by an earlier panel).
merge_focals    <- TRUE
merge_dist      <- 3
max_panel_terms <- 25

plot_w        <- 14        # figure width (inches; ~20% narrower than 16)
row_height_in <- 0.22        # vertical inches per GO-term row (height scaling)
panel_pad_in  <- 1.7         # fixed inches per panel (title + top gene labels + bottom FE axis)
tree_scale    <- 0.8         # tree width: x-units per nesting level (lower = narrower)
tree_gap      <- 0.2         # gap (x-units) between tree tips and the FE bars
bar_w_units   <- 6           # x data-units spanning the FE bar axis (0 .. fe_axis_max)
bar_gap       <- 0.8         # gap (x-units) between the FE bars and the gene tiles

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

# --- gene-category colours ---------------------------------------------------
cat_levels <- c("network gene", "dog seed", "human seed",
                "network & dog seed", "network & human seed")
# Logical scheme: dog = red, human = blue (fixed, used elsewhere). A seed gene
# that is ALSO in the network gets a darker shade of its own seed colour, so the
# combined categories read as "emphasized" versions of the seed. Network-only is
# a distinct neutral hue (purple). No grey.
cat_colors <- c("network gene"          = "#c39bd3",   # network only (light purple)
                "dog seed"              = "#e1706c",   # dog seed only (red)
                "human seed"            = "#5c80bd",   # human seed only (blue)
                "network & dog seed"    = "#9e3b34",   # dog + network (dark red)
                "network & human seed"  = "#36568c")   # human + network (dark blue)

## ---------------------------------------------------------------------------
## LOAD SEED GENES (frozen CSV)
## ---------------------------------------------------------------------------
# If the sheet is shared with "anyone with the link", no login is required:
# If the sheet is PRIVATE, comment the line above and authenticate instead, e.g.

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
  message("Note: ", length(both_seed), " gene(s) flagged as BOTH dog and human seed; ",
          "they will be coloured as the dog category (priority).")
message("Seed genes: ", length(dog_seed_genes), " dog, ", length(human_seed_genes), " human")

## ---------------------------------------------------------------------------
## GO-TERM -> human gene SYMBOLS (authoritative GOA annotations, DIRECT only)
##   Direct gene -> GO annotations are read from the GO Consortium human GAF.
##   A gene belongs to a term iff it is annotated to THAT EXACT term in the GAF:
##   annotations are NOT propagated to descendant terms (and `regulates` cross-
##   links are not followed). This is stricter than org.Hs.egGO2ALLEGS. Used to
##   decide whether a seed gene belongs to a plotted GO term even when it was
##   not part of the enrichment query (network) intersection.
## ---------------------------------------------------------------------------

if (!file.exists(goa_gaf)) {
  message("Downloading authoritative GO annotations:\n  ", goa_url)
  utils::download.file(goa_url, goa_gaf, mode = "wb", quiet = TRUE)
}
message("Reading GO annotations from ", goa_gaf, " ...")
# GAF 2.2: tab-delimited, '!' comment/header lines. Columns of interest:
#   3 = DB Object Symbol, 4 = Qualifier, 5 = GO ID, 9 = Aspect (P/F/C).
gaf <- readr::read_tsv(goa_gaf, comment = "!", col_names = FALSE, quote = "",
                       col_types = readr::cols(.default = "c"), progress = FALSE)
gaf_direct <- gaf %>%
  dplyr::filter(X9 %in% c("P", "F", "C"),
                !grepl("NOT", X4, fixed = TRUE)) %>%   # drop negated annotations
  dplyr::transmute(symbol = norm_sym(X3), go = X5) %>%
  dplyr::distinct()
# term -> directly-annotated symbols (no propagation; see genes_in_term below)
direct2genes <- split(gaf_direct$symbol, gaf_direct$go)
message("GOA: ", nrow(gaf_direct), " direct annotations across ",
        length(direct2genes), " GO terms")

genes_in_term <- local({
  cache <- new.env(parent = emptyenv())
  function(go_id) {
    if (is.null(go_id) || is.na(go_id) || !grepl("^GO:", go_id))
      return(character(0))                  # non-GO sources (KEGG/REAC/...) -> none
    if (!is.null(cache[[go_id]])) return(cache[[go_id]])
    # DIRECT annotations only: a gene belongs to a term iff it is annotated to
    # that exact term in the GAF (no propagation to descendants / regulates).
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
    # fold enrichment = (k/n) / (K/N)
    .fe     = (.isize / .data[[col_qsize]]) /
      (.data[[col_termsize]] / .data[[col_effdom]])
  ) %>%
  # drop terms with fewer than min_overlap overlapping genes (applies to every
  # plotted row: focals and their ancestor/sibling relatives alike)
  dplyr::filter(.isize >= min_overlap)
message(nrow(go_ds), " rows remain after requiring >= ", min_overlap,
        " overlapping genes")

## ---------------------------------------------------------------------------
## RESOLVE COMMUNITIES
##   If no list is supplied, use the n_biggest clusters. Cluster size is taken
##   as the query_size (number of genes submitted to enrichment for that
##   community), which is constant within a community.
## ---------------------------------------------------------------------------
if (is.null(communities) || length(communities) == 0) {
  comm_size <- go_ds %>%
    dplyr::group_by(.comm) %>%
    dplyr::summarise(size = max(.data[[col_qsize]], na.rm = TRUE),
                     .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(size))
  communities <- head(comm_size$.comm, n_biggest)
  message("No communities supplied; using the ", length(communities),
          " biggest cluster(s) by query_size: ",
          paste0(head(comm_size$.comm, n_biggest),
                 " (n=", head(comm_size$size, n_biggest), ")",
                 collapse = ", "))
}
stopifnot(length(communities) > 0)

## ---------------------------------------------------------------------------
## GENE x COMMUNITY OVERLAP MATRIX
##   Writes ONE tab-delimited file (across all communities): one row per gene
##   that appears in any community, one column per community, TRUE/FALSE for
##   membership. A gene counts as "in" a community if it appears in any
##   eligible significant GO set for that community (the same `go_ds` pool the
##   tile plots are built from). Rows are ordered by how many communities a
##   gene is shared across (descending), then alphabetically.
## ---------------------------------------------------------------------------
genes_by_comm <- setNames(
  lapply(communities, function(cm)
    sort(unique(unlist(go_ds$.genes[go_ds$.comm == cm])))),
  communities
)

empty_comm <- communities[lengths(genes_by_comm) == 0]
if (length(empty_comm))
  message("Note: no eligible genes for ", paste(empty_comm, collapse = ", "),
          " (all-FALSE column[s]).")

all_overlap_genes <- sort(unique(unlist(genes_by_comm)))
if (length(all_overlap_genes) == 0) {
  message("No genes found across the requested communities; overlap file not written.")
} else {
  overlap_mat <- tibble::tibble(gene = all_overlap_genes)
  for (cm in communities)
    overlap_mat[[cm]] <- overlap_mat$gene %in% genes_by_comm[[cm]]
  
  # order: most-shared genes first, then alphabetical
  ord         <- order(-rowSums(as.matrix(overlap_mat[communities])), overlap_mat$gene)
  overlap_mat <- overlap_mat[ord, , drop = FALSE]
  
  overlap_tsv <- paste0(plot_ds, "_community_gene_overlap.tsv")
  readr::write_tsv(overlap_mat, overlap_tsv)
  message("Wrote ", overlap_tsv, " (", nrow(overlap_mat), " genes x ",
          length(communities), " communities)")
}

## ---------------------------------------------------------------------------
## SHARED HELPERS
## ---------------------------------------------------------------------------
get_safe <- function(id, map) tryCatch(AnnotationDbi::get(id, map),
                                       error = function(e) character(0))
san      <- function(x) str_replace_all(x, ":", "_")

# wrap a long label onto two lines, breaking at the space closest to the middle
wrap_mid <- function(s, max_chars = 80) {
  vapply(s, function(x) {
    if (is.na(x) || nchar(x) <= max_chars) return(x)
    sp <- gregexpr(" ", x, fixed = TRUE)[[1]]
    if (length(sp) == 1 && sp[1] == -1) return(x)   # no space to break on
    brk <- sp[which.min(abs(sp - nchar(x) / 2))]    # space nearest the middle
    paste0(substr(x, 1, brk - 1), "\n", substr(x, brk + 1, nchar(x)))
  }, character(1), USE.NAMES = FALSE)
}

# Classify one (gene, GO term) tile into one of the 5 categories.
# `net_syms` = uppercased network (query-intersection) symbols for the term.
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

# rectangular cladogram coordinates for an ape phylo (no extra packages):
# returns node, parent, x (depth from root), y (tip order), isTip, label
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
##   Returns list(focal, onto, sel) or NULL if the focal is absent.
## ---------------------------------------------------------------------------
compute_sel <- function(fg) {
  
  if (!(fg %in% sig_terms)) {
    warning("Focal term ", fg, " is not in the eligible pool - skipping.")
    return(NULL)
  }
  
  onto <- str_extract(sub$.source[sub$.go == fg][1], "BP|MF|CC")
  if (is.na(onto))
    onto <- tryCatch(as.character(AnnotationDbi::Ontology(GO.db::GOTERM[[fg]])),
                     error = function(e) NA_character_)
  if (is.na(onto) || !onto %in% c("BP", "MF", "CC")) {
    warning("Cannot resolve GO ontology for ", fg, " - skipping.")
    return(NULL)
  }
  ANCESTOR <- switch(onto, BP = GO.db::GOBPANCESTOR, MF = GO.db::GOMFANCESTOR, CC = GO.db::GOCCANCESTOR)
  PARENTS  <- switch(onto, BP = GO.db::GOBPPARENTS,  MF = GO.db::GOMFPARENTS,  CC = GO.db::GOCCPARENTS)
  CHILDREN <- switch(onto, BP = GO.db::GOBPCHILDREN, MF = GO.db::GOMFCHILDREN, CC = GO.db::GOCCCHILDREN)
  
  sel <- fg
  if (include_ancestors) {
    sel <- union(sel, intersect(setdiff(get_safe(fg, ANCESTOR), "all"), sig_terms))
  } else {
    sel <- union(sel, intersect(get_safe(fg, PARENTS), sig_terms))
  }
  if (include_siblings) {
    par  <- get_safe(fg, PARENTS)
    sibs <- unique(unlist(lapply(par, get_safe, map = CHILDREN)))
    sel  <- union(sel, intersect(sibs, sig_terms))
  }
  if (include_children) {
    sel <- union(sel, intersect(get_safe(fg, CHILDREN), sig_terms))
  }
  sel <- unique(c(fg, intersect(sel, sig_terms)))
  message("[", fg, "] ", onto, ": ", length(sel), " GO terms")
  list(focal = fg, onto = onto, sel = sel)
}

## ---------------------------------------------------------------------------
## STEP 2: build the GO-nesting tree, then one tile panel per focal
## ---------------------------------------------------------------------------
build_tree <- function(sel, onto) {
  
  ANCESTOR <- switch(onto, BP = GO.db::GOBPANCESTOR, MF = GO.db::GOMFANCESTOR, CC = GO.db::GOCCANCESTOR)
  
  ## reduce the GO DAG to a tree: attach each term to its nearest in-set ancestor;
  ## give every term a self-leaf so internal ancestors also get a row
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
  
  to_newick <- function(node) {
    kids <- children_map[[node]]
    if (length(kids) == 0) return(san(node))                      # pure leaf
    parts <- c(san(node), vapply(kids, to_newick, character(1)))  # self-leaf + kids
    paste0("(", paste(parts, collapse = ","), ")")                # unlabelled internal
  }
  if (length(sel) == 1) {
    nwk <- paste0("(", san(sel), ");")
  } else if (length(roots) == 1) {
    nwk <- paste0(to_newick(roots), ";")
  } else {
    nwk <- paste0("(", paste(vapply(roots, to_newick, character(1)),
                             collapse = ","), ");")
  }
  ape::read.tree(text = nwk)
}

## ---------------------------------------------------------------------------
## Build one panel from a prebuilt tree. `tile_x0` (the x of the first gene
## column) is shared across all panels in a PDF so the columns line up.
## ---------------------------------------------------------------------------
build_panel <- function(sel, title, tree, tile_x0, pad_to = NULL, fe_limits = NULL) {
  
  lay     <- tree_layout(tree)
  lay$x   <- lay$x * tree_scale          # compress the tree horizontally
  xmax    <- max(lay$x)
  
  ## plot data
  sel_sub    <- sub %>% dplyr::filter(.go %in% sel)
  net_by_go  <- setNames(sel_sub$.genes, sel_sub$.go)   # network (query) genes per term
  fe_by_go   <- setNames(sel_sub$.fe,    sel_sub$.go)
  name_by_go <- setNames(sel_sub$.name,  sel_sub$.go)
  pval_by_go <- setNames(sel_sub$.pval,  sel_sub$.go)
  
  # uppercased network symbols per term (for case-insensitive classification)
  netU_by_go <- lapply(net_by_go, norm_sym)
  
  # seed genes directly annotated to each term that are NOT already network
  # genes -> these are the extra "dog seed"/"human seed" tiles.
  seedextra_by_go <- setNames(lapply(sel, function(g) {
    ann      <- genes_in_term(g)                       # uppercased human symbols
    extra    <- setdiff(intersect(ann, c(dog_seed_genes, human_seed_genes)),
                        netU_by_go[[g]])
    extra
  }), sel)
  
  # displayed genes per term = network genes (original case) + extra seed genes.
  # Network labels keep their original case; seed-only labels are uppercase.
  disp_by_go <- setNames(lapply(sel, function(g) {
    unique(c(net_by_go[[g]], seedextra_by_go[[g]]))
  }), sel)
  
  all_genes  <- unique(unlist(disp_by_go))
  if (length(all_genes) == 0) all_genes <- character(0)
  gene_freq  <- table(factor(unlist(disp_by_go), levels = all_genes))
  # column order: network genes first (left), then seed-only (dog/human, never
  # a network gene in this panel). Within each block: frequency desc, then A-Z.
  net_all    <- unique(unlist(netU_by_go))            # all network symbols (uppercase)
  is_net_col <- norm_sym(all_genes) %in% net_all
  gene_order <- all_genes[order(!is_net_col,
                                -as.integer(gene_freq[all_genes]),
                                all_genes)]
  
  # one blank spacer column separates the network block from the seed-only block
  n_net_cols <- sum(is_net_col)
  n_genes    <- length(gene_order)
  col_pos    <- seq_len(n_genes)
  if (n_net_cols > 0 && n_net_cols < n_genes)
    col_pos[(n_net_cols + 1):n_genes] <- col_pos[(n_net_cols + 1):n_genes] + 1
  col_x      <- setNames(col_pos, gene_order)            # gene -> slot index (with gap)
  n_slots    <- if (n_genes) max(col_pos) else 0         # real slots incl. spacer
  
  # pad the right with blank columns so every panel has >= pad_to gene columns
  pad_to     <- pad_to %||% n_genes
  n_pad      <- max(0, pad_to - n_genes)
  n_col      <- n_slots + n_pad                          # total x slots
  
  # tip y-position per GO term (tips are labelled with sanitised ids)
  tip_y     <- setNames(lay$y[lay$isTip], lay$label[lay$isTip])
  
  # tile data: one row per (GO term, displayed gene), categorised + carrying FE
  tiledf <- do.call(rbind, lapply(sel, function(g) {
    gs <- intersect(disp_by_go[[g]], gene_order)
    if (!length(gs)) return(NULL)
    cats <- vapply(norm_sym(gs), classify_gene, character(1),
                   net_syms = netU_by_go[[g]])
    data.frame(x = tile_x0 + col_x[gs] - 0.5, y = tip_y[san(g)],
               fe = fe_by_go[[g]], go = g, gene = gs, category = cats,
               stringsAsFactors = FALSE)
  }))
  if (!is.null(tiledf))
    tiledf$category <- factor(tiledf$category, levels = cat_levels)
  
  # background grid: a very light grey outline for every cell (incl. empties)
  bg_grid <- expand.grid(go = sel, gene = gene_order,
                         KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  if (nrow(bg_grid)) {
    bg_grid$x <- tile_x0 + col_x[bg_grid$gene] - 0.5
    bg_grid$y <- unname(tip_y[san(bg_grid$go)])
    bg_grid   <- bg_grid[!is.na(bg_grid$y), , drop = FALSE]
  }
  
  # tree edges (rectangular): horizontal child segments + vertical parent spans
  ed   <- lay[!is.na(lay$parent), ]
  px   <- lay$x[match(ed$parent, lay$node)]
  py   <- lay$y[match(ed$parent, lay$node)]
  hseg <- data.frame(x = px,        xend = ed$x, y = ed$y, yend = ed$y)
  vseg <- do.call(rbind, lapply(unique(ed$parent), function(pp) {
    ys <- ed$y[ed$parent == pp]
    data.frame(x = lay$x[lay$node == pp], xend = lay$x[lay$node == pp],
               y = min(ys), yend = max(ys))
  }))
  
  # y-axis labels: star input focals, plus p value and fold enrichment (1 dp)
  star      <- ifelse(sel %in% focal_go, "★ ", "")
  disp      <- setNames(
    paste0(star, sel, "  ", name_by_go[sel],
           "  (p=", formatC(pval_by_go[sel], format = "e", digits = 1), ")"),
    san(sel))
  disp      <- setNames(wrap_mid(disp, 80), names(disp))
  tip_lab   <- lay[lay$isTip, ]
  tip_lab   <- tip_lab[order(tip_lab$y), ]
  
  # x-axis labels: gene names at their slots; blank for the spacer + padding
  x_breaks  <- tile_x0 + seq_len(n_col) - 0.5
  x_labels  <- rep("", n_col)
  x_labels[col_pos] <- gene_order
  
  # darker grey box, slightly enlarged, around the non-network (seed-only) gene
  # columns at the right end of the panel.
  seed_cols <- if (n_net_cols < n_genes) gene_order[(n_net_cols + 1):n_genes] else character(0)
  box_df <- NULL
  if (length(seed_cols) > 0 && length(tip_y) > 0) {
    s_slots <- col_x[seed_cols]
    box_df <- data.frame(
      xmin = tile_x0 + min(s_slots) - 1 - 0.18,  # left edge of first seed column
      xmax = tile_x0 + max(s_slots) + 0.18,      # right edge of last seed column
      ymin = min(tip_y) - 0.65,
      ymax = max(tip_y) + 0.65
    )
  }
  
  # ---- fold-change bar plot, to the LEFT of the gene tiles ----------------
  # Bars share a common FE axis (fe_axis_max) across every panel in the PDF so
  # bar lengths are directly comparable. FE=0 baseline sits at bar_x0 (just
  # right of the tree) and bars grow rightward toward the tiles. The FE tick
  # labels and "Fold enrichment" title are drawn UNDERNEATH the plot.
  fe_axis_max <- if (!is.null(fe_limits) && is.finite(max(fe_limits, na.rm = TRUE)))
    max(fe_limits, na.rm = TRUE) else suppressWarnings(max(fe_by_go[sel], na.rm = TRUE))
  if (!is.finite(fe_axis_max) || fe_axis_max <= 0) fe_axis_max <- 1
  bar_x0      <- tile_x0 - bar_gap - bar_w_units      # FE=0 baseline (left edge)
  
  bardf <- data.frame(y  = unname(tip_y[san(sel)]),
                      fe = unname(fe_by_go[sel]),
                      stringsAsFactors = FALSE)
  bardf <- bardf[!is.na(bardf$y) & !is.na(bardf$fe), , drop = FALSE]
  bardf$xmin <- bar_x0
  bardf$xmax <- bar_x0 + (bardf$fe / fe_axis_max) * bar_w_units
  
  # shared FE axis: gridlines through the bar block, ticks + title underneath
  fe_breaks <- pretty(c(0, fe_axis_max), n = 4)
  fe_breaks <- fe_breaks[fe_breaks >= 0 & fe_breaks <= fe_axis_max + 1e-9]
  bar_y_top <- if (length(tip_y)) max(tip_y) + 0.5 else 1
  bar_y_bot <- if (length(tip_y)) min(tip_y) - 0.5 else 0
  axis_y    <- bar_y_bot - 0.55                       # FE tick-label row (below)
  title_y   <- bar_y_bot - 1.15                       # "Fold enrichment" caption
  gx        <- bar_x0 + (fe_breaks / fe_axis_max) * bar_w_units
  grid_df   <- data.frame(x = gx, y = bar_y_bot, yend = bar_y_top)
  tick_df   <- data.frame(x = gx, y = axis_y, lab = formatC(fe_breaks, format = "g"))
  
  # per-tile long table for the CSV export
  long <- if (is.null(tiledf)) {
    tibble::tibble(panel = character(0), go = character(0), gene = character(0),
                   fe = numeric(0), category = character(0))
  } else {
    tibble::tibble(panel = title, go = tiledf$go, gene = tiledf$gene,
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
    # fold-change bars (shared scale): gridlines + dark-grey bars; FE axis below
    geom_segment(data = grid_df, aes(x = x, y = y, xend = x, yend = yend),
                 colour = "grey88", linewidth = 0.3, inherit.aes = FALSE) +
    geom_rect(data = bardf, aes(xmin = xmin, xmax = xmax,
                                ymin = y - 0.2, ymax = y + 0.2),
              fill = "grey30", colour = NA, inherit.aes = FALSE) +
    geom_text(data = tick_df, aes(x = x, y = y, label = lab),
              size = 1.9, colour = "grey40", vjust = 1, inherit.aes = FALSE) +
    annotate("text", x = bar_x0 + bar_w_units / 2, y = title_y,
             label = "Fold enrichment", size = 2.4, colour = "grey30", vjust = 1) +
    scale_fill_manual(name = "Gene type", values = cat_colors,
                      breaks = cat_levels, drop = FALSE,
                      na.value = "white") +
    scale_x_continuous(position = "top", breaks = x_breaks, labels = x_labels,
                       expand = expansion(mult = 0)) +
    scale_y_continuous(breaks = tip_lab$y, labels = disp[tip_lab$label],
                       expand = expansion(add = c(1.5, 0.5))) +
    coord_cartesian(xlim = c(-0.3, tile_x0 + n_col + 0.3), clip = "off") +
    labs(x = NULL, y = NULL, title = title) +
    theme_minimal(base_size = 9) +
    theme(
      axis.text.x     = element_text(angle = 90, hjust = 0, vjust = 0.5, size = 6),
      axis.text.y     = element_text(size = 7),
      panel.grid      = element_blank(),
      legend.position = "right",
      plot.title      = element_text(size = 10, face = "bold")
    )
  
  # outline the seed-only gene block (drawn on top so the border stays visible)
  if (!is.null(box_df))
    gg <- gg + geom_rect(data = box_df,
                         aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                         inherit.aes = FALSE, fill = NA,
                         colour = "grey62", linewidth = 0.6)
  
  list(gg = gg, n = length(sel), long = long)
}

## ---------------------------------------------------------------------------
## STEP 3: for each community, auto-pick the interesting focal terms (every term
##         with p < pvalue_max and intersection > min_intersection). Each focal
##         seeds a neighbourhood (focal + significant ancestors/siblings). When
##         merge_focals is TRUE, focals within merge_dist GO-hierarchy steps of
##         each other are grouped (closest first) into shared panels -- one tree
##         per group -- but no panel is allowed to exceed max_panel_terms rows,
##         so a large related cluster is split into a few bounded panels rather
##         than one giant one; no focal is dropped. When FALSE, the original
##         greedy behaviour is used (one panel per focal, dropping focals already
##         covered). One PDF per community.
## ---------------------------------------------------------------------------
out_pdfs <- character(0)

for (cm in communities) {
  
  # this community's eligible pool (rows that may appear: focals + relatives)
  sub <- go_ds %>%
    dplyr::filter(.comm == cm) %>%
    dplyr::group_by(.go) %>% dplyr::slice(1) %>% dplyr::ungroup()
  if (nrow(sub) == 0) {
    message("== ", cm, ": no eligible rows - skipping"); next
  }
  sig_terms <- sub$.go
  fe_of     <- setNames(sub$.fe, sub$.go)
  
  # auto-select interesting focal terms: every term passing the criteria,
  # ordered by fold enrichment (no cap on how many)
  focal_go <- sub %>%
    dplyr::filter(.pval < pvalue_max, .isize > min_intersection) %>%
    dplyr::arrange(dplyr::desc(.fe)) %>%
    dplyr::pull(.go)
  if (length(focal_go) == 0) {
    message("== ", cm, ": no interesting focal terms (p < ", pvalue_max,
            ", intersection > ", min_intersection, ") - skipping"); next
  }
  message("== ", cm, ": ", length(focal_go), " candidate focal terms ==")
  
  # one neighbourhood (focal + significant ancestors/siblings) per focal
  sels <- Filter(Negate(is.null), lapply(focal_go, compute_sel))
  if (length(sels) == 0) { message("  no resolvable focals for ", cm); next }
  nfoc <- length(sels)
  
  if (merge_focals) {
    # ---- merge nearby focals into fewer panels, bounded by a size cap --------
    sel_list <- lapply(sels, `[[`, "sel")
    onto_v   <- vapply(sels, `[[`, character(1), "onto")
    foc_v    <- vapply(sels, `[[`, character(1), "focal")
    
    # ancestor -> min #steps (up to merge_dist) above a focal, in its ontology
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
    
    # candidate links: same ontology, GO distance (steps from each focal up to
    # the nearest shared ancestor, summed) <= merge_dist; sorted closest-first
    cand <- list()
    if (nfoc > 1)
      for (i in 1:(nfoc - 1)) for (j in (i + 1):nfoc) {
        if (onto_v[i] != onto_v[j]) next
        common <- intersect(names(depl[[i]]), names(depl[[j]]))
        if (!length(common)) next
        d <- min(depl[[i]][common] + depl[[j]][common])
        if (is.finite(d) && d <= merge_dist) cand[[length(cand) + 1]] <- c(i, j, d)
      }
    
    # agglomerative union, closest first, but never exceed max_panel_terms rows
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
          parent[ri]    <- rj
          grp_sel[[rj]] <- merged
        }
      }
    }
    groups <- split(seq_len(nfoc), vapply(seq_len(nfoc), froot, integer(1)))
    
    keep <- lapply(unname(groups), function(ix) {
      foci <- foc_v[ix]
      sel  <- intersect(unique(unlist(sel_list[ix])), sig_terms)
      list(sel   = sel,
           onto  = onto_v[ix[1]],
           foci  = foci,
           genes = unique(unlist(sub$.genes[sub$.go %in% sel])),
           maxfe = suppressWarnings(max(fe_of[foci], na.rm = TRUE)),
           label = foci[which.max(fe_of[foci])])
    })
    keep <- keep[order(-vapply(keep, `[[`, numeric(1), "maxfe"))]
    message("  ", nfoc, " focal term(s) -> ", length(keep),
            " panel(s) (merge_dist=", merge_dist, ", cap=", max_panel_terms, " rows)")
  } else {
    # ---- original greedy: one panel per focal, drop already-covered focals ----
    sels    <- sels[order(-fe_of[vapply(sels, `[[`, character(1), "focal")])]
    covered <- character(0); keep <- list()
    for (s in sels) {
      if (length(intersect(s$sel, covered)) > 0) {
        message("  skip ", s$focal, " (FE=", signif(fe_of[s$focal], 3),
                ") - folded into a previous plot"); next
      }
      s$foci  <- s$focal
      s$label <- s$focal
      s$genes <- unique(unlist(sub$.genes[sub$.go %in% s$sel]))
      keep[[length(keep) + 1]] <- s
      covered <- union(covered, s$sel)
    }
  }
  
  # pad-to width: network genes plus any extra seed genes annotated to the
  # panel's terms, so columns line up and seed-only tiles always have a slot.
  panel_width <- function(s) {
    net  <- unique(unlist(sub$.genes[sub$.go %in% s$sel]))
    netU <- norm_sym(net)
    extra <- setdiff(intersect(unique(unlist(lapply(s$sel, genes_in_term))),
                               c(dog_seed_genes, human_seed_genes)), netU)
    length(unique(c(net, extra)))
  }
  max_genes <- max(vapply(keep, panel_width, numeric(1)))
  fe_limits <- range(fe_of[unlist(lapply(keep, `[[`, "sel"))], na.rm = TRUE)
  
  # build trees once, then reserve a common tree width so the first gene column
  # starts at the same x in every panel (-> columns line up across panels)
  trees   <- lapply(keep, function(s) build_tree(s$sel, s$onto))
  depths  <- vapply(trees, function(t) max(tree_layout(t)$x), numeric(1))
  tile_x0 <- max(depths) * tree_scale + tree_gap + bar_w_units + bar_gap
  message("  building ", length(keep), " panel(s); padding to ", max_genes,
          " gene columns; tiles start at x=", tile_x0)
  
  panels <- Map(function(s, tr) {
    title <- paste0(plot_ds, " ", cm, ": ", s$label,
                    if (length(s$foci) > 1) " & related focal terms" else " & overlapping terms")
    build_panel(s$sel, title, tree = tr, tile_x0 = tile_x0,
                pad_to = max_genes, fe_limits = fe_limits)
  }, keep, trees)
  
  n_rows     <- vapply(panels, `[[`, numeric(1), "n")
  heights_in <- n_rows * row_height_in + panel_pad_in   # per-panel height (inches)
  
  # Draw the legend only ONCE for the whole PDF, built from a dummy plot that
  # contains EVERY category, so all values render with their colour (absent
  # categories otherwise show as blank/white keys in a panel-derived legend).
  legend_src <- ggplot(
    data.frame(cat = factor(cat_levels, levels = cat_levels),
               x = 1, y = seq_along(cat_levels)),
    aes(x, y, fill = cat)) +
    geom_tile() +
    scale_fill_manual(name = "Gene type", values = cat_colors,
                      breaks = cat_levels, drop = FALSE) +
    theme_minimal(base_size = 9) +
    theme(legend.position = "right", legend.box.margin = margin(0, 0, 0, 6))
  legend <- cowplot::get_legend(legend_src)
  
  ggs    <- lapply(panels, function(p) p$gg + theme(legend.position = "none"))
  body   <- cowplot::plot_grid(
    plotlist    = ggs,
    ncol        = 1,
    rel_heights = heights_in,
    align       = "v",         # align panels vertically...
    axis        = "lr"         # ...by their left and right edges
  )
  # Legend pinned to the TOP of the right-hand column. Its cell gets a fixed
  # height in inches (sized to the number of keys) rather than a fraction of the
  # figure, so short figures still reserve enough room for the full key instead
  # of clipping it. The blank cell below absorbs the remaining height, and the
  # figure height is floored so it is never shorter than the legend itself.
  legend_in <- 0.22 * length(cat_levels) + 0.8        # natural legend height (in)
  fig_h     <- max(sum(heights_in), legend_in + 0.3)  # never shorter than legend
  right_col <- cowplot::plot_grid(legend, NULL, ncol = 1,
                                  rel_heights = c(legend_in, fig_h - legend_in))
  combined  <- cowplot::plot_grid(body, right_col, ncol = 2, rel_widths = c(1, 0.16))
  
  out_stem <- paste0(plot_ds, "_", cm, "_p", -log10(pvalue_max))
  out_pdf  <- paste0(out_stem, "_nesting_tileplot.pdf")
  out_csv  <- paste0(out_stem, "_nesting_tileplot_data.csv")
  ggsave(out_pdf, plot = combined, width = plot_w,
         height = fig_h, limitsize = FALSE)
  readr::write_csv(
    dplyr::bind_rows(lapply(panels, `[[`, "long")) %>%
      dplyr::transmute(community = cm, panel, go, gene,
                       fold_enrichment = fe, category),
    out_csv
  )
  out_pdfs <- c(out_pdfs, out_pdf)
  message("  wrote ", out_pdf, " (", length(panels), " panel(s), ",
          sum(n_rows), " rows)")
}

message("Done. ", length(out_pdfs), " PDF(s): ", paste(out_pdfs, collapse = ", "))