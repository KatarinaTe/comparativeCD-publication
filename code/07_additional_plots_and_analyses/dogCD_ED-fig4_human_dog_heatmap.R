library(Seurat)
#BiocManager::install("fgsea")
library(fgsea)
library(tidyverse)
#BiocManager::install("UCell")
library(UCell)
library(ComplexHeatmap)
library(circlize)

#### Comparison of expression of OCD/dogCD network genes between human and dog cells ########

# Load in the integrated human-dog data set

dh_5regions <- readRDS("./human_dog_5regions.rds")

cell_types <- unique(dh_5regions$supercluster_term)


# This creates a matrix showing cell counts for Species (rows) by Supercluster (columns)
table(dh_5regions$species, dh_5regions$supercluster_term, useNA = "ifany")

# Calculate average expression per species

avg_expr <- AverageExpression(dh_5regions, 
                              group.by = c("supercluster_term", "species"), 
                              features = dogCD_OCD_network, 
                              assays = "RNA")$RNA

# To get a quantitative measure of how well the overall expression profile of these 165 genes is conserved within each cell type, calculate the Spearman correlation.

correlation_results <- data.frame(Cell_Type = character(), Correlation = numeric(), P_Value = numeric())

for (ct in cell_types) {
  human_col <- paste0(ct, "_human")
  dog_col <- paste0(ct, "_dog")
  
  # Check if both species have this cell type
  if (human_col %in% colnames(avg_expr) & dog_col %in% colnames(avg_expr)) {
    test <- cor.test(avg_expr[, human_col], avg_expr[, dog_col], method = "spearman")
    correlation_results <- rbind(correlation_results, 
                                 data.frame(Cell_Type = ct, 
                                            Correlation = test$estimate, 
                                            P_Value = test$p.value))
  }
}

print(correlation_results)

# Note: High positive correlations indicate the 165 genes behave very similarly in that specific cell type across both species.

# 1. Separate the expression matrices
human_cols <- grep("_human$", colnames(avg_expr), value = TRUE)
dog_cols <- grep("_dog$", colnames(avg_expr), value = TRUE)

mat_human <- as.matrix(avg_expr[, human_cols])
mat_dog <- as.matrix(avg_expr[, dog_cols])

# Ensure column orders match (e.g., strip the species suffix and sort)
colnames(mat_human) <- gsub("_human$", "", colnames(mat_human))
colnames(mat_dog) <- gsub("_dog$", "", colnames(mat_dog))

common_genes <- intersect(rownames(mat_human), rownames(mat_dog))
mat_human <- mat_human[common_genes, ]
mat_dog <- mat_dog[common_genes, ]

# Define neuron cell types

neurons <- c("UL-tel","DL-tel","DL-cor-6b","DL-NP","Th-ex","Cer-ex","Cer-i","Mid-i","CGE-inter","MGE-inter","LAMP5")

# 1. Identify cell types that exist in BOTH human and dog matrices
# (Assuming mat_human and mat_dog from the previous steps)
common_long_names <- intersect(colnames(mat_human), colnames(mat_dog))

# Subset the matrices to only include these common cells
mat_human_common <- mat_human[, common_long_names]
mat_dog_common <- mat_dog[, common_long_names]

# 2. Rename the columns from long names to abbreviations
# unname() ensures we just pass the string, not the dictionary key
colnames(mat_human_common) <- unname(abbrev_map[colnames(mat_human_common)])
colnames(mat_dog_common) <- unname(abbrev_map[colnames(mat_dog_common)])

# 3. Filter your cell_abbrevs to only include cells that actually exist in the common set
# This prevents errors if a mapped cell type is missing from the data
valid_order <- intersect(cell_abbrevs, colnames(mat_human_common))

# 4. Reorder both matrices to match your exact specified order
mat_human_ordered <- mat_human_common[, valid_order]
mat_dog_ordered <- mat_dog_common[, valid_order]

# Calculate Correlation
cor_matrix_ordered <- cor(mat_human_ordered, mat_dog_ordered, method = "spearman")

# Create a grouping factor to separate Neurons from Glia/Other
# We check if the row/col name is in our 'neurons' list
cell_groupings <- ifelse(rownames(cor_matrix_ordered) %in% neurons, "Neurons", "Glia/Other")

# Convert to a factor and set levels so "Neurons" always prints first (top/left)
split_factor <- factor(cell_groupings, levels = c("Neurons", "Glia/Other"))

# Plot the Final Heatmap
Heatmap(cor_matrix_ordered,
        name = "Spearman\nCorrelation",
        col = col_fun_cor, # (Uses the color ramp from previous step)
        
        # CRITICAL: Turn off clustering to enforce your custom order
        cluster_rows = FALSE,    
        cluster_columns = FALSE, 
        
        # Split the heatmap into Neurons and Glia
        row_split = split_factor,
        column_split = split_factor,
        
        # Formatting titles
        row_title = "Human",
        column_title = "Dog",
        row_title_rot = 0, # Keeps "Neurons" and "Glia" text horizontal
        
        # Add values inside the boxes (optional)
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", cor_matrix_ordered[i, j]), x, y, 
                    gp = gpar(fontsize = 8, col = ifelse(cor_matrix_ordered[i, j] > 0.7, "white", "black")))
        })

