#!/usr/bin/env Rscript
# cross_species_expression.R — Fig 4a-b, dogCD/OCD network gene-set enrichment against the
# Siletti et al. 2023 human brain snRNA-seq atlas, by cell type and by brain region.
#
# Source: network_gene_expression_enrichment_plots_code_updated.r (collaborator script, left
# untouched at its original location for provenance) — literal reproduction, only the
# readRDS() path replaced with a CLI arg, and explicit saveRDS()/ggsave() calls added for the
# 2 derived matrices and 2 plots the original built in-memory but never wrote to disk (needed
# so Nextflow has real output files to capture — no other logic changed).
#
# Dropped 2026-09-22: library(ComplexHeatmap)/library(circlize), loaded in the original but never
# actually called anywhere in this script (the Spearman-correlation heatmap they'd be for is a
# separate, not-yet-converted analysis — see Supplementary Note 8). Confirmed via package
# changelogs that this is the only container consumer affected by the Seurat/fgsea version
# question raised in REPRODUCIBILITY_AUDIT.md's 2026-09-22 tool-version update; dropping these
# two removes ComplexHeatmap/circlize's version from that question entirely, since they're unused.
suppressPackageStartupMessages({
  library(Seurat)
  library(fgsea)
  library(UCell)
  library(dplyr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
siletti_rds <- args[1]

# Network Gene set
dogCD_OCD_network <- c("CDH12", "CDH18", "CDH9", "LINGO1", "LRRC4C", "PAK3", "PLXNA4", "SEMA6C", "CDH10", "CDH8", "EPHA4", "GAP43", "NTRK2", "RBFOX1", "SEMA6D", "DAAM2",
                       "NTRK3", "PCDH19", "RHOB", "SEMA6B", "SRGAP3", "CHL1", "CNTN6", "L1CAM", "NFASC", "CNTN1", "CNTNAP1", "LSAMP", "NCAM2", "NRCAM", "DLGAP1", "LRRTM1",
                       "CACNA1H", "CACNB2", "HCN2", "KCNB1", "KCND3", "KCNV1", "BSN", "CACNG3", "GRIA2", "GRIA3", "GRIA4", "GRIK1", "GRIK2", "GRIK3", "GRIN2A", "GRIN2B", "GRM3",
                       "GRM5", "GRM7", "KCNC1", "LRRTM2", "SCN2A", "STXBP1", "SYT4", "ADGRB3", "PDE1B", "PDE4D", "SEPTIN3", "ARIH1", "ARIH2", "UBA7", "UBE2L6", "AJAP1", "CADM2",
                       "CNTNAP4", "LRFN2", "LRFN3", "LRFN5", "NYNRIN", "PCSK1N", "SLC4A10", "ATP1A2", "DGKI", "ELK3", "GALK2", "GAS7", "GATAD2A", "OMG", "SCN2B", "AKAP6", "CELF5",
                       "CRTC1", "DACH1", "DAOA", "DCLK1", "DPP10", "ELAVL2", "ELAVL3", "GABRA1", "GABRA2", "GABRA5", "GABRA6", "GABRB2", "GABRB3", "GABRG1", "KCNIP1", "LRRN1", "NLGN1",
                       "NOVA1", "NRGN", "OPCML", "PPFIA2", "ADGRA1", "ADGRD1", "ADGRF3", "ADGRG2", "ADGRG7", "ANKIB1", "BIK", "C20orf144", "C2CD3", "CA9", "CCHCR1", "ENC1", "EVI2B", "F2RL1",
                       "FHOD1", "GPR137", "HLA-DMA", "IGSF11", "IQCC", "IQSEC3", "IVNS1ABP", "KANSL2", "KATNB1", "KATNIP", "KBTBD12", "KLHDC8B", "KLHL1", "KLHL22", "KLHL24", "KLHL26", "KLHL31",
                       "KLHL33", "KLHL35", "KLHL38", "KLHL40", "KLHL41", "KLHL5", "KLHL6", "KLHL7", "MTSS2", "NECAB3", "NES", "NPHP4", "NR2C2", "NRM", "P4HTM", "PAWR", "PLAAT4", "PROSER3",
                       "PWWP2A", "RASGRF1", "RBL1", "RNASEH2A", "SAMD9L", "SIGIRR", "SPMIP2", "SPSB4", "TEX9", "UNC80", "ZNF396", "ZNF839")

# Process the Siletti human dataset to calculate per-cell-type average expression
human_data <- readRDS(siletti_rds)

# Create named list of genes to check
features_list <- list(genes_list = dogCD_OCD_network)

# 2. Calculate UCell scores and append them to Seurat metadata
# (Creates a new column in metadata called "TargetGenes_UCell")
human_data <- AddModuleScore_UCell(human_data, features = features_list)

# Define neurons and glial cells in metadata
# 1. Define the categories
neurons <- c("UL-tel", "DL-tel", "DL-cor-6b", "DL-NP", "Am-ex", "Th-ex",
             "LRL", "URL", "Cer-i", "Mid-i", "CGE-inter", "MGE-inter",
             "LAMP5", "EMS", "MS", "H-CA1-3", "H-CA4", "H-DG", "MB", "Splatter")

glial <- c("Astro", "Bergmann", "Choroid", "Ependymal", "Fibro", "Micro",
           "Oligo", "OP", "COP", "Vascular")

# 2. Create the new metadata column
human_data$cell_class <- case_when(
  human_data$cell_type_abbrev %in% neurons ~ "Neuron",
  human_data$cell_type_abbrev %in% glial ~ "Glial",
  human_data$cell_type_abbrev == "Misc" ~ "Other",
  TRUE ~ "Unknown" # Catch-all safely handles any unexpected or NA values
)

# 3. Verify the new metadata
table(human_data$cell_class, human_data$cell_type_abbrev)

# 3. Extract metadata
human_metadata <- human_data@meta.data
colnames(human_metadata)

# Retain only neurons and glial cells
clean_metadata <- human_metadata %>%
  filter(cell_class %in% c("Neuron", "Glial"))

# Verify the filtering worked
table(clean_metadata$cell_class, useNA = "always")

# 4. Perform a statistical test
# Direct comparison via Wilcoxon Rank-Sum test
human_CD_wilcox_test <- wilcox.test(genes_list_UCell ~ cell_class, data = clean_metadata)
print(human_CD_wilcox_test)

# Look at the summary stats
clean_metadata %>%
  group_by(cell_class) %>%
  summarize(
    mean_expression = mean(genes_list_UCell, na.rm = TRUE),
    median_expression = median(genes_list_UCell, na.rm = TRUE),
    cell_count = n() # Always good practice to see how many cells are in each group
  )

# 2. Get the average expression per cell type
# This returns a matrix where rows are genes and columns are your cell types
human_avg_by_celltype <- AverageExpression(
  human_data,
  group.by = "cell_type_abbrev",
  slot = "data" # Use "data" for log-normalized values or "counts" for raw
)$RNA

# 3. Get the global average expression across ALL cells
# We create a temporary dummy column in the metadata to group everything together
human_data$total_pool <- "all"
human_avg_global <- AverageExpression(
  human_data,
  group.by = "total_pool",
  slot = "data" # Use "data" for log-normalized values or "counts" for raw
)$RNA


# 4. Calculate the specificity ratio
# We divide the cell-type matrix row-wise by the single global average vector
human_cell_specificity_scores <- human_avg_by_celltype / human_avg_global[, 1]

# 5. Clean up edge cases
# If a gene is completely unexpressed anywhere, 0/0 will result in NaN.
human_cell_specificity_scores[is.nan(human_cell_specificity_scores)] <- 0

# Check your new specificity matrix
head(human_cell_specificity_scores)

# Save as RDS object
# Use standard saveRDS, but explicit compression can optimize space
saveRDS(human_cell_specificity_scores, file = "human_cell_specificity_scores.rds", compress = "xz")

# Read it back
human_cell_specificity_scores <- readRDS("human_cell_specificity_scores.rds")


# Now look at the full set of network genes in the human data

# fgsea expects a named list of gene sets, even if you only have one
human_network_pathways <- list(network = dogCD_OCD_network)

# 2. Create an empty list to store results for each cell type
human_network_gsea_results_list <- list()

# 3. Loop through each cell type (column) in your specificity matrix
for (celltype in colnames(human_cell_specificity_scores)) {

  # Extract the scores for the current cell type
  scores <- human_cell_specificity_scores[, celltype]

  # fgsea requires a named numeric vector sorted in descending order
  ranked_genes <- sort(scores, decreasing = TRUE)

  # Run fgsea
  # Note: 'minSize = 2' ensures it won't crash if your small subset
  # only has a couple of genes matching your dataset.
  fgsea_res <- fgsea(
    pathways = human_network_pathways,
    stats = ranked_genes,
    minSize = 2,
    maxSize = Inf
  )

  # Add a column so we know which cell type these results belong to
  fgsea_res$cell_type <- celltype

  # Save to our list
  human_network_gsea_results_list[[celltype]] <- fgsea_res
}

# 4. Combine all cell type results into one clean data frame
human_network_gsea_results <- do.call(rbind, human_network_gsea_results_list)

# 5. View the results sorted by the strongest enrichment
human_network_gsea_results <- human_network_gsea_results[order(human_network_gsea_results$pval), ]
head(human_network_gsea_results)

# Define cells into neurons and glial cells

human_network_gsea_results <- human_network_gsea_results %>%
  mutate(cell_class = case_when(
    # List all neuron abbreviations here
    cell_type %in% c("UL-tel","DL-tel","DL-cor-6b","DL-NP","Am-ex","Th-ex","LRL","URL","Cer-i","Mid-i","CGE-inter","MGE-inter","LAMP5","EMS","MS","H-CA1-3","H-CA4","H-DG","MB","Splatter") ~ "Neurons",

    # List all your glial abbreviations here
    cell_type %in% c("Astro","Bergmann","Choroid","Ependymal","Fibro","Micro","Oligo","OP","COP","Vascular") ~ "Glial cells",

    # Fallback option for anything else (e.g., Endothelial, Pericytes)
    TRUE ~ "Other"
  ))

# plot the p-values

human_network_gsea_results$cell_class <- factor(human_network_gsea_results$cell_class, levels = c("Neurons", "Glial cells", "Other"), ordered=TRUE)
#scale_fill_manual(values = c("orange","brown", "grey")) +

human_network_gsea_results_plot <- ggplot(data=human_network_gsea_results, aes(x=reorder(cell_type, log10(padj)), y = -log10(padj), fill = cell_class)) +
  geom_col() +
  geom_hline(yintercept = -log10(0.05), colour="tomato", lty=2) +
  geom_hline(yintercept = -log10(0.1), colour="purple", lty=2) +
  scale_fill_manual(values = c("orange","brown", "grey")) +
  theme_minimal() +
  labs(
    title = "",
    x = "Cell Type",
    y = "-log10(adj.P)"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 9),
    axis.text.y = element_text(size = 7), # Small text for 85 genes
    panel.grid = element_blank()
  )

# Nextflow output capture — the original script built this plot object but never saved it.
ggsave("human_network_gsea_results_plot.pdf", human_network_gsea_results_plot, width = 10, height = 6)


# And the same again but per brain region, and just for neurons

# First, create a subset that only contains neurons

human_neurons <- subset(human_data[,human_data$cell_type_abbrev %in% c("UL-tel","DL-tel","DL-cor-6b","DL-NP","Am-ex","Th-ex","LRL","URL","Cer-i","Mid-i","CGE-inter","MGE-inter","LAMP5","EMS","MS","H-CA1-3","H-CA4","H-DG","MB","Splatter")])

# 2. Get the average expression per cell type
# This returns a matrix where rows are genes and columns are your cell types
human_neurons_avg_by_tissue <- AverageExpression(
  human_neurons,
  group.by = "ROIGroupFine",
  slot = "data" # Use "data" for log-normalized values or "counts" for raw
)$RNA

# and across all regions

human_neurons_avg_global <- AverageExpression(
  human_neurons,
  group.by = "total_pool",
  slot = "data" # Use "data" for log-normalized values or "counts" for raw
)$RNA

# 4. Calculate the specificity ratio
# We divide the tissue-type matrix row-wise by the single global average vector
human_tissue_specificity_scores <- human_neurons_avg_by_tissue / human_neurons_avg_global[, 1]

# 5. Clean up edge cases
# If a gene is completely unexpressed anywhere, 0/0 will result in NaN.
human_tissue_specificity_scores[is.nan(human_tissue_specificity_scores)] <- 0

# Check your new specificity matrix
head(human_tissue_specificity_scores)

# Nextflow output capture — the original script built this matrix but never saved it, unlike
# its cell-type counterpart above.
saveRDS(human_tissue_specificity_scores, file = "human_tissue_specificity_scores.rds", compress = "xz")

### Now for the gene set enrichment analysis - for each tissue type, are genes in our GWAS gene set over-represented in the most tissue-type specific genes?

human_network_pathways <- list(network = dogCD_OCD_network)

# 2. Create an empty list to store results for each cell type
human_network_per_region_gsea_results_list <- list()

# 3. Loop through each tissue type (column) in your specificity matrix
for (tissuetype in colnames(human_tissue_specificity_scores)) {

  # Extract the scores for the current tissue type
  scores <- human_tissue_specificity_scores[, tissuetype]

  # fgsea requires a named numeric vector sorted in descending order
  ranked_genes <- sort(scores, decreasing = TRUE)

  # Run fgsea
  # Note: 'minSize = 2' ensures it won't crash if your small subset
  # only has a couple of genes matching your dataset.
  fgsea_res <- fgsea(
    pathways = human_network_pathways,
    stats = ranked_genes,
    minSize = 2,
    maxSize = Inf
  )

  # Add a column so we know which tissue type these results belong to
  fgsea_res$tissue_type <- tissuetype

  # Save to our list
  human_network_per_region_gsea_results_list[[tissuetype]] <- fgsea_res
}

# 4. Combine all tissue type results into one clean data frame
human_network_per_region_gsea_results <- do.call(rbind, human_network_per_region_gsea_results_list)


# 5. View the results sorted by the strongest enrichment
human_network_per_region_gsea_results <- human_network_per_region_gsea_results[order(human_network_per_region_gsea_results$pval), ]
head(human_network_per_region_gsea_results)


# plot the p-values

human_network_per_region_gsea_results_plot <- ggplot(data=human_network_per_region_gsea_results, aes(x=reorder(tissue_type, log10(padj)), y = -log10(padj))) +
  geom_col(fill="grey") +
  geom_hline(yintercept = -log10(0.05), colour="tomato", lty=2) +
  geom_hline(yintercept = -log10(0.1), colour="purple", lty=2) +
  theme_minimal() +
  labs(
    title = "",
    x = "Human brain region",
    y = "-log10 p value"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 9),
    axis.text.y = element_text(size = 7), # Small text for 85 genes
    panel.grid = element_blank()
  )

# Nextflow output capture — the original script built this plot object but never saved it.
ggsave("human_network_per_region_gsea_results_plot.pdf", human_network_per_region_gsea_results_plot, width = 10, height = 6)
