#!/usr/bin/env Rscript
# cross_species_transcriptomic_heatmap.R -- Extended Data Fig. 3: heatmap of pairwise Spearman
# rank correlation coefficients between human and dog pseudo-bulk expression profiles for the
# 165 OCD/dogCD network genes, across matched cell types (neurons vs. glia/other).
#
# Source: ../dogCD_ED-fig4_human_dog_heatmap.R (collaborator script, left untouched at its
# original location for provenance) -- literal reproduction, with these changes:
#   - readRDS() path and the network gene list replaced with CLI args (the gene list was
#     hardcoded a second time in the collaborator's script with no visible source; it's
#     identical to the 165-gene data/05_network/ocd_ccd_systemsmap_genes.txt already used by
#     the sibling 3_cross_species_expression stage, so read from there instead of a third
#     hardcoded copy).
#   - abbrev_map/cell_abbrevs/col_fun_cor -- referenced in the original but never defined
#     anywhere in it or found anywhere else in this repo -- filled in below. See each one's own
#     comment for exactly how it was resolved, and this stage's README for the full table.
#   - `cell_types` narrowed to the 20 supercluster_term values actually plotted in the published
#     figure (the original used all 26) -- 2026-09-30, per KatarinaTe: "plot what is plotted in
#     ED fig 3, nothing more, and the same with spearman only for those to sum up in the end."
#     One-line change, same loop/structure as the original otherwise: AverageExpression() still
#     runs the same way over the full object, and the correlation loop is still the original's
#     own cor.test-per-cell-type pattern over avg_expr's columns -- just narrower scope, not a
#     restructure.
#   - explicit saveRDS()/pdf() calls added for the outputs the original built in-memory but
#     never wrote to disk (needed so Nextflow has real output files to capture).
#   - AverageExpression() given slot = "data" explicitly (log-normalized values), matching the
#     sibling 3_cross_species_expression stage's own explicit choice -- the original script left
#     this as Seurat's default, which is version-sensitive; making it explicit avoids that.
suppressPackageStartupMessages({
  library(Seurat)
  library(ComplexHeatmap)
  library(circlize)
})

args <- commandArgs(trailingOnly = TRUE)
human_dog_rds <- args[1]
seed_gene_list <- args[2]
out_dir <- args[3]

dogCD_OCD_network <- readLines(seed_gene_list)

# ---- abbrev_map: long supercluster_term -> short axis label used in the published heatmap.
# Resolved 2026-09-30 by reading the real embedded heatmap image (Extended Data Fig. 3, from
# NATURE_Extended_Data_Figures_260929.docx's image2.png) axis-by-axis and matching each label
# back to the real supercluster_term values pulled from this .rds. 9 of the 20 (the glia/
# vascular/other block) confirmed directly by KatarinaTe. The 11 neuron labels match this
# repo's OWN already-converted, manuscript-verified sibling script's `neurons` vector
# (3_cross_species_expression/bin/cross_species_expression.R) exactly, including "URL" as its
# abbreviation for the developmental lineage producing excitatory cerebellar granule neurons in
# the Siletti atlas's terminology -- independent corroboration for the one label ("Cer-ex" here
# vs "URL" there -- same underlying supercluster_term, "Upper rhombic lip") that wasn't
# self-evident from the image alone. See this stage's README for the full one-by-one table.
abbrev_map <- c(
  "Upper-layer intratelencephalic" = "UL-tel",
  "Deep-layer intratelencephalic" = "DL-tel",
  "Deep-layer corticothalamic and 6b" = "DL-cor-6b",
  "Deep-layer near-projecting" = "DL-NP",
  "Thalamic excitatory" = "Th-ex",
  "Upper rhombic lip" = "Cer-ex",
  "Cerebellar inhibitory" = "Cer-i",
  "Midbrain-derived inhibitory" = "Mid-i",
  "CGE interneuron" = "CGE-inter",
  "MGE interneuron" = "MGE-inter",
  "LAMP5-LHX6 and Chandelier" = "LAMP5",
  "Astrocyte" = "Astro",
  "Bergmann glia" = "Bergmann",
  "Ependymal" = "Ependymal",
  "Fibroblast" = "Fibro",
  "Microglia" = "Micro",
  "Oligodendrocyte" = "Oligo",
  "Oligodendrocyte precursor" = "OP",
  "Committed oligodendrocyte precursor" = "COP",
  "Vascular" = "Vascular"
)

# Axis display order, read directly off the published heatmap image (neurons block, then
# glia/vascular/other block, left-to-right / top-to-bottom).
cell_abbrevs <- c("UL-tel", "DL-tel", "DL-cor-6b", "DL-NP", "Th-ex", "Cer-ex", "Cer-i", "Mid-i",
                   "CGE-inter", "MGE-inter", "LAMP5",
                   "Astro", "Bergmann", "Ependymal", "Fibro", "Micro", "Oligo", "OP", "COP", "Vascular")

neurons <- c("UL-tel", "DL-tel", "DL-cor-6b", "DL-NP", "Th-ex", "Cer-ex", "Cer-i", "Mid-i",
             "CGE-inter", "MGE-inter", "LAMP5")

glia <- setdiff(cell_abbrevs, neurons)

# ---- col_fun_cor: also referenced but never defined in the original. Not recoverable from the
# data (it's a plotting choice, not biological data). Approximated here with the standard
# 'viridis' palette's own well-known anchor hex colours, in their normal (non-reversed) order --
# low correlations dark (purple/blue), high correlations light (yellow) -- matching the published
# heatmap. (Briefly reversed 2026-09-30, switched back 2026-10-01 per KatarinaTe: the original
# figure has yellow = higher correlation, blue = lower.) Rather than pulling in a new package
# dependency for one colour ramp. Revisit if the collaborator provides their exact function --
# doesn't affect any of the correlation numbers below, only the heatmap's colour rendering.
col_fun_cor <- colorRamp2(
  seq(0, 1, length.out = 5),
  c("#440154FF", "#3B528BFF", "#21908CFF", "#5DC863FF", "#FDE725FF")
)

# Load the integrated human-dog-mouse Seurat object, mouse cells already dropped by the
# collaborator ("I generated it from the integrated human-dog-mouse single cell Seurat object I
# made, after dropping all mouse cells from the dataset").
dh_5regions <- readRDS(human_dog_rds)

# This creates a matrix showing cell counts for Species (rows) by Supercluster (columns)
print(table(dh_5regions@meta.data$species, dh_5regions@meta.data$supercluster_term, useNA = "ifany"))

# Only the 20 supercluster_term values actually plotted in the published figure -- narrowed from
# the original's `unique(dh_5regions$supercluster_term)` (all 26). AverageExpression() computes
# each group's mean independently of which other groups are present, so this narrowing changes
# nothing about the numbers for these 20 -- it only stops the correlation loop below from also
# reporting the 6 unplotted types.
cell_types <- intersect(unique(dh_5regions@meta.data$supercluster_term), names(abbrev_map))

# Calculate average expression per species
avg_expr <- AverageExpression(dh_5regions,
                              group.by = c("supercluster_term", "species"),
                              features = dogCD_OCD_network,
                              assays = "RNA",
                              slot = "data")$RNA

# To get a quantitative measure of how well the overall expression profile of these 165 genes is
# conserved within each cell type, calculate the Spearman correlation.
correlation_results <- data.frame(Cell_Type = character(), Correlation = numeric(), P_Value = numeric())
for (ct in cell_types) {
  human_col <- paste0(ct, "_human")
  dog_col <- paste0(ct, "_dog")
  if (human_col %in% colnames(avg_expr) & dog_col %in% colnames(avg_expr)) {
    test <- cor.test(avg_expr[, human_col], avg_expr[, dog_col], method = "spearman")
    correlation_results <- rbind(correlation_results,
                                 data.frame(Cell_Type = ct, Abbrev = unname(abbrev_map[ct]),
                                            Correlation = unname(test$estimate), P_Value = test$p.value))
  }
}
write.csv(correlation_results, file.path(out_dir, "human_dog_correlation_by_celltype.csv"), row.names = FALSE)
print(correlation_results)

# Note: High positive correlations indicate the 165 genes behave very similarly in that specific
# cell type across both species.

# Separate the expression matrices
human_cols <- grep("_human$", colnames(avg_expr), value = TRUE)
dog_cols <- grep("_dog$", colnames(avg_expr), value = TRUE)

mat_human <- as.matrix(avg_expr[, human_cols])
mat_dog <- as.matrix(avg_expr[, dog_cols])

colnames(mat_human) <- gsub("_human$", "", colnames(mat_human))
colnames(mat_dog) <- gsub("_dog$", "", colnames(mat_dog))

common_genes <- intersect(rownames(mat_human), rownames(mat_dog))
mat_human <- mat_human[common_genes, ]
mat_dog <- mat_dog[common_genes, ]

# Restrict to cell types present in both species, rename to abbreviations, order to match the
# published figure
common_long_names <- intersect(colnames(mat_human), colnames(mat_dog))
mat_human_common <- mat_human[, common_long_names]
mat_dog_common <- mat_dog[, common_long_names]

colnames(mat_human_common) <- unname(abbrev_map[colnames(mat_human_common)])
colnames(mat_dog_common) <- unname(abbrev_map[colnames(mat_dog_common)])

valid_order <- intersect(cell_abbrevs, colnames(mat_human_common))
mat_human_ordered <- mat_human_common[, valid_order]
mat_dog_ordered <- mat_dog_common[, valid_order]

# Calculate Correlation (full pairwise matrix -- needed for the heatmap's off-diagonal cells too)
cor_matrix_ordered <- cor(mat_human_ordered, mat_dog_ordered, method = "spearman")
saveRDS(cor_matrix_ordered, file.path(out_dir, "human_dog_correlation_matrix.rds"))

cell_groupings <- ifelse(rownames(cor_matrix_ordered) %in% neurons, "Neurons", "Glia/Other")
split_factor <- factor(cell_groupings, levels = c("Neurons", "Glia/Other"))

pdf(file.path(out_dir, "human_dog_heatmap.pdf"), width = 10, height = 9)
draw(Heatmap(cor_matrix_ordered,
        name = "Spearman\nCorrelation",
        col = col_fun_cor,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        row_split = split_factor,
        column_split = split_factor,
        row_title = "Human",
        column_title = "Dog",
        row_title_rot = 0,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", cor_matrix_ordered[i, j]), x, y,
                    gp = gpar(fontsize = 8, col = ifelse(cor_matrix_ordered[i, j] > 0.7, "white", "black")))
        }))
dev.off()

# Matched-cell-type Spearman r for the neuron block specifically -- the number to check against
# the manuscript's "Cross-species correlations are highest in matching neuronal cell types
# (Spearman r = 0.80-0.89)".
cat("\nMatched-cell-type Spearman r, neuron block (published range: 0.80-0.89):\n")
print(correlation_results[correlation_results$Abbrev %in% neurons, ])

# Same summary for the glia/vascular/other block, for comparison -- no published range given for
# this block specifically, just useful to see whether it's meaningfully lower than the neuron
# block, as the heatmap's own colour scale would suggest.
cat("\nMatched-cell-type Spearman r, glia/vascular/other block:\n")
print(correlation_results[correlation_results$Abbrev %in% glia, ])
