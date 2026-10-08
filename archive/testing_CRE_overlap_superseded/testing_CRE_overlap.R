library(tidyverse)

# Counts of brain and other tissue CREs from EPIC data
N_brain <- 126410 
N_other <- 238469 

# Analysis of overlap between clump +- 100Kb regions and EPIC dog CREs - do we see significant enrichment of active brain elements
# in these regions?


# 1. Set up the contingency table for testing significance
# O_brain  = Overlap with GWAS (Brain)
O_brain <- 408   
    
# O_other = Overlap with GWAS (Other)
O_other <- 655       
  
# Calculate the "No Overlap" groups
no_O_brain <- N_brain - O_brain
no_O_other <- N_other - O_other

# 2. Construct the 2x2 Contingency Table
# byrow = TRUE ensures the data fills row-by-row mapping to our dimnames
contingency_table <- matrix(
  c(O_brain, no_O_brain,
    O_other, no_O_other),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(
    Tissue = c("Brain", "Other"),
    GWAS_Overlap = c("Yes", "No")
  )
)

# 3. Run the Fisher's Exact Test
# We use alternative="greater" because you want to know if Brain > Other. 
# If you just wanted to know if they were *different*, you'd use "two.sided" (the default).
fisher_results <- fisher.test(contingency_table, alternative = "greater")

# 4. View the results
print(fisher_results)

# Answer: Yes! but not a huge effect: odds ratio = 1.17564, p = 0.0059

## Compare this to the SIZE regions

size_O_brain <- 325   
size_O_other <- 736       

# Calculate the "No Overlap" groups
size_no_O_brain <- N_brain - size_O_brain
size_no_O_other <- N_other - size_O_other

# 2. Construct the 2x2 Contingency Table
# byrow = TRUE ensures the data fills row-by-row mapping to our dimnames
size_contingency_table <- matrix(
  c(size_O_brain, size_no_O_brain,
    size_O_other, size_no_O_other),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(
    Tissue = c("Brain", "Other"),
    GWAS_Overlap = c("Yes", "No")
  )
)

# 3. Run the Fisher's Exact Test
# We use alternative="greater" because you want to know if Brain > Other. 
size_fisher_results <- fisher.test(size_contingency_table, alternative = "greater")

# 4. View the results
print(size_fisher_results)

# Answer: Not significant! p = 0.9975, OR = 0.83

## Ok, now test for the GWS plus suggestive SNPs clumps

# First, for dogCD
# 1. Set up the contingency table for testing significance
dogCD_GWS_sugg_O_brain <- 698  
dogCD_GWS_sugg_O_other <- 1300

# Calculate the "No Overlap" groups
dogCD_GWS_sugg_no_O_brain <- N_brain - dogCD_GWS_sugg_O_brain
dogCD_GWS_sugg_no_O_other <- N_other - dogCD_GWS_sugg_O_other

# 2. Construct the 2x2 Contingency Table
# byrow = TRUE ensures the data fills row-by-row mapping to our dimnames
dogCD_GWS_sugg_contingency_table <- matrix(
  c(dogCD_GWS_sugg_O_brain, dogCD_GWS_sugg_no_O_brain,
    dogCD_GWS_sugg_O_other, dogCD_GWS_sugg_no_O_other),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(
    Tissue = c("Brain", "Other"),
    GWAS_Overlap = c("Yes", "No")
  )
)

# 3. Run the Fisher's Exact Test
# We use alternative="greater" because you want to know if Brain > Other. 
# If you just wanted to know if they were *different*, you'd use "two.sided" (the default).
dogCD_GWS_sugg_fisher_results <- fisher.test(dogCD_GWS_sugg_contingency_table, alternative = "greater")

# 4. View the results
print(dogCD_GWS_sugg_fisher_results)

# Answer:

## Compare this to the SIZE regions

size_GWS_sugg_O_brain <- 430   
size_GWS_sugg_O_other <- 996       


# Calculate the "No Overlap" groups
size_GWS_sugg_no_O_brain <- N_brain - size_GWS_sugg_O_brain
size_GWS_sugg_no_O_other <- N_other - size_GWS_sugg_O_other

# 2. Construct the 2x2 Contingency Table
# byrow = TRUE ensures the data fills row-by-row mapping to our dimnames
size_GWS_sugg_contingency_table <- matrix(
  c(size_GWS_sugg_O_brain, size_GWS_sugg_no_O_brain,
    size_GWS_sugg_O_other, size_GWS_sugg_no_O_other),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(
    Tissue = c("Brain", "Other"),
    GWAS_Overlap = c("Yes", "No")
  )
)

# 3. Run the Fisher's Exact Test
# We use alternative="greater" because you want to know if Brain > Other. 
# If you just wanted to know if they were *different*, you'd use "two.sided" (the default).
size_GWS_sugg_fisher_results <- fisher.test(size_GWS_sugg_contingency_table, alternative = "greater")

# 4. View the results
print(size_GWS_sugg_fisher_results)

# Answer: Not significant! p = 0.9999, OR = 0.81

######################


## Ok, let's redo this but rather than counting elements (which are highly variable in size), we count the base pair overlap, i.e. how many bps of our
# GWAS regions is covered by CREs
# Need to count the total bps of our brain and other regions CREs:

brain_cre_bp <- 247096600
other_cre_bp <- 393644600

# Now let's look at the DogCD GWS clumps +/- 100KB

dogCD_GWS_O_bp_brain <- 799607  
dogCD_GWS_O_bp_other <- 1109740      

# Calculate the "No Overlap" groups
dogCD_GWS_no_O_bp_brain <- brain_cre_bp - dogCD_GWS_O_bp_brain
dogCD_GWS_no_O_bp_other <- other_cre_bp - dogCD_GWS_O_bp_other

# 2. Construct the 2x2 Contingency Table
# byrow = TRUE ensures the data fills row-by-row mapping to our dimnames
dogCD_GWS_bp_contingency_table <- matrix(
  c(dogCD_GWS_O_bp_brain, dogCD_GWS_no_O_bp_brain,
    dogCD_GWS_O_bp_other, dogCD_GWS_no_O_bp_other),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(
    Tissue = c("Brain", "Other"),
    GWAS_Overlap = c("Yes", "No")
  )
)

# 3. Run the Fisher's Exact Test
# We use alternative="greater" because you want to know if Brain > Other. 
# If you just wanted to know if they were *different*, you'd use "two.sided" (the default).
dogCD_GWS_bp_fisher_results <- fisher.test(dogCD_GWS_bp_contingency_table, alternative = "greater")

# 4. View the results
print(dogCD_GWS_bp_fisher_results)

# Answer: p < 2.2x10-16, OR = 1.15

# Now let's look at the SIZE GWS clumps +/- 100KB

SIZE_GWS_O_bp_brain <- 710315
SIZE_GWS_O_bp_other <- 1330316

# Calculate the "No Overlap" groups
SIZE_GWS_no_O_bp_brain <- brain_cre_bp - SIZE_GWS_O_bp_brain
SIZE_GWS_no_O_bp_other <- other_cre_bp - SIZE_GWS_O_bp_other

# 2. Construct the 2x2 Contingency Table
# byrow = TRUE ensures the data fills row-by-row mapping to our dimnames
SIZE_GWS_bp_contingency_table <- matrix(
  c(SIZE_GWS_O_bp_brain, SIZE_GWS_no_O_bp_brain,
    SIZE_GWS_O_bp_other, SIZE_GWS_no_O_bp_other),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(
    Tissue = c("Brain", "Other"),
    GWAS_Overlap = c("Yes", "No")
  )
)

# 3. Run the Fisher's Exact Test
# We use alternative="greater" because you want to know if Brain > Other. 
# If you just wanted to know if they were *different*, you'd use "two.sided" (the default).
SIZE_GWS_bp_fisher_results <- fisher.test(SIZE_GWS_bp_contingency_table, alternative = "greater")

# 4. View the results
print(SIZE_GWS_bp_fisher_results)

# Answer: p = 1, OR = 0.85

#######################
#Now adding in UU regions:::
## Ok, let's redo this but rather than counting elements (which are highly variable in size), we count the base pair overlap, i.e. how many bps of our
# GWAS regions is covered by CREs
# Need to count the total bps of our brain and other regions CREs:

brain_cre_bp <- 247096600

# Counts when including UU regions
brain_cre_bp <- 460303000

other_cre_bp <- 393644600

# Now let's look at the DogCD GWS clumps +/- 100KB

dogCD_GWS_O_bp_brain <- 799607

# If including UU regions
dogCD_GWS_O_bp_brain <- 1472848

dogCD_GWS_O_bp_other <- 1109740      

# Calculate the "No Overlap" groups
dogCD_GWS_no_O_bp_brain <- brain_cre_bp - dogCD_GWS_O_bp_brain
dogCD_GWS_no_O_bp_other <- other_cre_bp - dogCD_GWS_O_bp_other

# 2. Construct the 2x2 Contingency Table
# byrow = TRUE ensures the data fills row-by-row mapping to our dimnames
dogCD_GWS_bp_contingency_table <- matrix(
  c(dogCD_GWS_O_bp_brain, dogCD_GWS_no_O_bp_brain,
    dogCD_GWS_O_bp_other, dogCD_GWS_no_O_bp_other),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(
    Tissue = c("Brain", "Other"),
    GWAS_Overlap = c("Yes", "No")
  )
)

# 3. Run the Fisher's Exact Test
# We use alternative="greater" because you want to know if Brain > Other. 
# If you just wanted to know if they were *different*, you'd use "two.sided" (the default).
dogCD_GWS_bp_fisher_results <- fisher.test(dogCD_GWS_bp_contingency_table, alternative = "greater")

# 4. View the results
print(dogCD_GWS_bp_fisher_results)

# Answer: p < 2.2x10-16, OR = 1.15

# Now let's look at the SIZE GWS clumps +/- 100KB

SIZE_GWS_O_bp_brain <- 710315

#With uu regions
SIZE_GWS_O_bp_brain <- 1272105

SIZE_GWS_O_bp_other <- 1330316

# Calculate the "No Overlap" groups
SIZE_GWS_no_O_bp_brain <- brain_cre_bp - SIZE_GWS_O_bp_brain
SIZE_GWS_no_O_bp_other <- other_cre_bp - SIZE_GWS_O_bp_other

# 2. Construct the 2x2 Contingency Table
# byrow = TRUE ensures the data fills row-by-row mapping to our dimnames
SIZE_GWS_bp_contingency_table <- matrix(
  c(SIZE_GWS_O_bp_brain, SIZE_GWS_no_O_bp_brain,
    SIZE_GWS_O_bp_other, SIZE_GWS_no_O_bp_other),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(
    Tissue = c("Brain", "Other"),
    GWAS_Overlap = c("Yes", "No")
  )
)

# 3. Run the Fisher's Exact Test
# We use alternative="greater" because you want to know if Brain > Other. 
# If you just wanted to know if they were *different*, you'd use "two.sided" (the default).
SIZE_GWS_bp_fisher_results <- fisher.test(SIZE_GWS_bp_contingency_table, alternative = "greater")

# 4. View the results
print(SIZE_GWS_bp_fisher_results)

# Answer: p = 1, OR = 0.85





## How about if we reduce the regions to just the 'clumps' regions without the +-100Kb

# 1. Set up the contingency table for testing significance
# O_brain  = Overlap with GWAS (Brain)
# N_brain  = Total Brain elements
O_brain <- 98   
N_brain <- 126410     

# O_other = Overlap with GWAS (Other)
# N_other = Total Other elements
O_other <- 173       
N_other <- 238469    

# Calculate the "No Overlap" groups
no_O_brain <- N_brain - O_brain
no_O_other <- N_other - O_other

# 2. Construct the 2x2 Contingency Table
# byrow = TRUE ensures the data fills row-by-row mapping to our dimnames
contingency_table <- matrix(
  c(O_brain, no_O_brain,
    O_other, no_O_other),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(
    Tissue = c("Brain", "Other"),
    GWAS_Overlap = c("Yes", "No")
  )
)

# 3. Run the Fisher's Exact Test
# We use alternative="greater" because you want to know if Brain > Other. 
# If you just wanted to know if they were *different*, you'd use "two.sided" (the default).
fisher_results <- fisher.test(contingency_table, alternative = "greater")

# 4. View the results
print(fisher_results)
# No longer significant, odds ration = 1.07, p = 0.32


## Ok, now we know there is significant enrichment (in the wider region set), lets test if there is significant enrichment within
# any of our specific brain regions in the UU CRE set using a chi-squared test

# 1. Bring in the overlap counts
UU_overlap_regions <- read.delim("~/Documents/GitHub/comparativeCD/data/05_network/cCRE/brain_region_overlaps.txt", sep =" ")

# Calculate No Overlap
UU_overlap_regions$No_overlap <- UU_overlap_regions$Total - UU_overlap_regions$Overlap

# Extract dynamically from the dataframe so you never have to type them out
regions <- UU_overlap_regions$Region

# Extract only the two columns we need and convert to a matrix
count_matrix <- as.matrix(UU_overlap_regions[, c("Overlap", "No_overlap")])

# Assign the brain regions as the row names of the matrix
rownames(count_matrix) <- regions

# 2. Now perform the Global Chi-Square Test
global_test <- chisq.test(count_matrix)
# N.S. p = 0.44

# What if we look at number of bp overlaps, rather than element counts...

# 1. Bring in the overlap counts
#/proj/caninebrain_2025/CCD_brain_CRE_intersect/CRE_overlap/UU_elements
UU_overlaps <- read.delim("~/Documents/GitHub/comparativeCD/data/05_network/cCRE/brain_region_bp_overlaps.txt", sep =" ")

# Calculate No Overlap
UU_overlaps$No_overlap_bp <- UU_overlaps$Total_bp - UU_overlaps$Overlap_bp

# Extract dynamically from the dataframe so you never have to type them out
regions <- UU_overlaps$Region

# Extract only the two columns we need and convert to a matrix
count_matrix <- as.matrix(UU_overlaps[, c("Overlap_bp", "No_overlap_bp")])

# Assign the brain regions as the row names of the matrix
rownames(count_matrix) <- regions
count_matrix_dogCD  <- count_matrix
# 2. Now perform the Global Chi-Square Test
global_test <- chisq.test(count_matrix)

# 3. "One vs. Rest" Fisher's Tests with FDR Correction
# We only proceed if the global test suggests a difference exists
if(global_test$p.value < 0.05) {
  cat("\nGlobal test significant! Running Post-Hoc 'One vs Rest' tests...\n")
  
  # Create empty vectors to store our raw p-values and odds ratios
  raw_p_values <- numeric(length(regions))
  odds_ratios <- numeric(length(regions))
  
  for(i in 1:length(regions)) {
    # 'One' is the current region - pulling directly from the dataframe columns!
    one_overlap <- UU_overlaps$Overlap_bp[i]
    one_no_overlap <- UU_overlaps$No_overlap_bp[i]
    
    # 'Rest' is the sum of all OTHER regions
    rest_overlap <- sum(UU_overlaps$Overlap_bp[-i])
    rest_no_overlap <- sum(UU_overlaps$No_overlap_bp[-i])
    
    # Build the 2x2 table for this specific test
    temp_table <- matrix(c(one_overlap, one_no_overlap, 
                           rest_overlap, rest_no_overlap), 
                         nrow = 2, byrow = TRUE)
    
    # Run one-sided Fisher's test (checking for enrichment: One > Rest)
    f_test <- fisher.test(temp_table, alternative = "greater")
    
    # Store results
    raw_p_values[i] <- f_test$p.value
    odds_ratios[i] <- f_test$estimate
  }
  
  # Apply Benjamini-Hochberg FDR correction to the p-values
  fdr_adjusted_p <- p.adjust(raw_p_values, method = "BH")
  
  # Compile everything into a neat results table
  results_df <- data.frame(
    Region = regions,
    Odds_Ratio = round(odds_ratios, 2),
    Raw_P = raw_p_values,
    FDR_Adjusted_P = fdr_adjusted_p
  )
  
  # Sort by most significant to least significant
  results_df <- results_df[order(results_df$FDR_Adjusted_P), ]
  
  cat("\n--- Final Enrichment Results ---\n")
  print(results_df)
  
} else {
  cat("\nGlobal test is not significant. No single region is driving the overlap more than others.\n")
}

# Plot the results
# 1. Calculate the absolute percentage of overlapping base pairs
UU_overlaps$Percent_Overlap <- (UU_overlaps$Overlap_bp / UU_overlaps$Total_bp) * 100

# 2. Reorder the 'Region' factor so the plot sorts from highest to lowest
# This makes it instantly visually clear which region is the top hit
UU_overlaps$Region <- factor(UU_overlaps$Region, 
                             levels = UU_overlaps$Region[order(UU_overlaps$Percent_Overlap, decreasing = TRUE)])

# 3. Create the bar plot
overlap_plot <- ggplot(UU_overlaps, aes(x = Region, y = Percent_Overlap, fill = Region)) +
  geom_bar(stat = "identity", color = "black", show.legend = FALSE) +
  
  # Add the exact percentage text on top of each bar
  geom_text(aes(label = sprintf("%.2f%%", Percent_Overlap)), 
            vjust = -0.5, size = 4.5) +
  
  # Expand the y-axis slightly so the text doesn't get cut off at the top
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  
  # Apply a clean, professional theme
  theme_minimal(base_size = 14) +
  labs(
    title = "Percentage of dogCD GWS Overlap by Brain Region",
    x = "Brain Region",
    y = "Percentage of Base Pairs Overlapping (%)"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    plot.title = element_text(face = "bold"),
    panel.grid.major.x = element_blank() # Removes vertical grid lines for a cleaner look
  )

# 4. Display the plot
print(overlap_plot)

###SIZE:##########
# 1. Bring in the overlap bp counts

UU_overlaps <- read.delim("~/Documents/GitHub/comparativeCD/data/05_network/cCRE/SIZE_brain_region_bp_overlaps.txt", sep ="\t")

# Calculate No Overlap
UU_overlaps$No_overlap_bp <- UU_overlaps$Total_bp - UU_overlaps$Overlap_bp

# Extract dynamically from the dataframe so you never have to type them out
regions <- UU_overlaps$Region

# Extract only the two columns we need and convert to a matrix
count_matrix <- as.matrix(UU_overlaps[, c("Overlap_bp", "No_overlap_bp")])
#

# Assign the brain regions as the row names of the matrix
rownames(count_matrix) <- regions
count_matrix_size  <- count_matrix

# 2. Now perform the Global Chi-Square Test
global_test <- chisq.test(count_matrix)

# 3. "One vs. Rest" Fisher's Tests with FDR Correction
# We only proceed if the global test suggests a difference exists
if(global_test$p.value < 0.05) {
  cat("\nGlobal test significant! Running Post-Hoc 'One vs Rest' tests...\n")
  
  # Create empty vectors to store our raw p-values and odds ratios
  raw_p_values <- numeric(length(regions))
  odds_ratios <- numeric(length(regions))
  
  for(i in 1:length(regions)) {
    # 'One' is the current region - pulling directly from the dataframe columns!
    one_overlap <- UU_overlaps$Overlap_bp[i]
    one_no_overlap <- UU_overlaps$No_overlap_bp[i]
    
    # 'Rest' is the sum of all OTHER regions
    rest_overlap <- sum(UU_overlaps$Overlap_bp[-i])
    rest_no_overlap <- sum(UU_overlaps$No_overlap_bp[-i])
    
    # Build the 2x2 table for this specific test
    temp_table <- matrix(c(one_overlap, one_no_overlap, 
                           rest_overlap, rest_no_overlap), 
                         nrow = 2, byrow = TRUE)
    
    # Run one-sided Fisher's test (checking for enrichment: One > Rest)
    f_test <- fisher.test(temp_table, alternative = "greater")
    
    # Store results
    raw_p_values[i] <- f_test$p.value
    odds_ratios[i] <- f_test$estimate
  }
  
  # Apply Benjamini-Hochberg FDR correction to the p-values
  fdr_adjusted_p <- p.adjust(raw_p_values, method = "BH")
  
  # Compile everything into a neat results table
  results_df <- data.frame(
    Region = regions,
    Odds_Ratio = round(odds_ratios, 2),
    Raw_P = raw_p_values,
    FDR_Adjusted_P = fdr_adjusted_p
  )
  
  # Sort by most significant to least significant
  results_df <- results_df[order(results_df$FDR_Adjusted_P), ]
  
  cat("\n--- Final Enrichment Results ---\n")
  print(results_df)
  
} else {
  cat("\nGlobal test is not significant. No single region is driving the overlap more than others.\n")
}



# Plot the results
# 1. Calculate the absolute percentage of overlapping base pairs
UU_overlaps$Percent_Overlap <- (UU_overlaps$Overlap_bp / UU_overlaps$Total_bp) * 100

# 2. Reorder the 'Region' factor so the plot sorts from highest to lowest
# This makes it instantly visually clear which region is the top hit
UU_overlaps$Region <- factor(UU_overlaps$Region, 
                             levels = UU_overlaps$Region[order(UU_overlaps$Percent_Overlap, decreasing = TRUE)])

# 3. Create the bar plot
overlap_plot <- ggplot(UU_overlaps, aes(x = Region, y = Percent_Overlap, fill = Region)) +
  geom_bar(stat = "identity", color = "black", show.legend = FALSE) +
  
  # Add the exact percentage text on top of each bar
  geom_text(aes(label = sprintf("%.2f%%", Percent_Overlap)), 
            vjust = -0.5, size = 4.5) +
  
  # Expand the y-axis slightly so the text doesn't get cut off at the top
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  
  # Apply a clean, professional theme
  theme_minimal(base_size = 14) +
  labs(
    title = "Percentage of SIZE GWAS Overlap by Brain Region",
    x = "Brain Region",
    y = "Percentage of Base Pairs Overlapping (%)"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    plot.title = element_text(face = "bold"),
    panel.grid.major.x = element_blank() # Removes vertical grid lines for a cleaner look
  )

# 4. Display the plot
print(overlap_plot)


####################################
#### compare dogCD vs size:
#total dogCD_gws_clump100kb region encompass: 6790724
#total size_gws_clump100kb region encompass: 6377108
doggenome = 2400000000 
dogCDbp = 6790724
Size_bp = 6377108

dogCDbp/doggenome * 100
Size_bp/doggenome * 100

count_matrix_size
count_matrix_dogCD

count_bp_dogCD <- as.data.frame(count_matrix_dogCD)
sum(count_bp_dogCD$Overlap_bp)
mean(count_bp_dogCD$Overlap_bp+count_bp_dogCD$No_overlap_bp)


##testing:
dogCDbp <- 6790724
Size_bp  <- 6377108

dogCD_overlap <- as.data.frame(count_matrix_dogCD)$Overlap_bp
size_overlap  <- as.data.frame(count_matrix_size)$Overlap_bp

regions <- rownames(count_matrix_dogCD)

results <- data.frame(
  region       = regions,
  rate_dogCD   = dogCD_overlap / dogCDbp,
  rate_size    = size_overlap  / Size_bp,
  OR           = NA,
  p            = NA
)

for (i in seq_along(regions)) {
  mat <- matrix(c(dogCD_overlap[i],  dogCDbp - dogCD_overlap[i],
                  size_overlap[i],   Size_bp  - size_overlap[i]),
                nrow = 2, byrow = TRUE,
                dimnames = list(c("dogCD","size"), c("overlap","no_overlap")))
  ft <- fisher.test(mat)
  results$OR[i] <- ft$estimate
  results$p[i]  <- ft$p.value
}

results$p_fdr <- p.adjust(results$p, method = "fdr")
print(results)


#####

dogCDbp <- 6790724
Size_bp  <- 6377108

dogCD_overlap <- as.data.frame(count_matrix_dogCD)$Overlap_bp
size_overlap  <- as.data.frame(count_matrix_size)$Overlap_bp

regions <- rownames(count_matrix_dogCD)

results <- data.frame(
  region       = regions,
  rate_dogCD   = dogCD_overlap / dogCDbp,
  rate_size    = size_overlap  / Size_bp,
  OR           = NA,
  CI_low       = NA,
  CI_high      = NA,
  p            = NA
)

for (i in seq_along(regions)) {
  mat <- matrix(c(dogCD_overlap[i],  dogCDbp - dogCD_overlap[i],
                  size_overlap[i],   Size_bp  - size_overlap[i]),
                nrow = 2, byrow = TRUE,
                dimnames = list(c("dogCD","size"), c("overlap","no_overlap")))
  ft <- fisher.test(mat)
  results$OR[i]      <- ft$estimate
  results$CI_low[i]  <- ft$conf.int[1]
  results$CI_high[i] <- ft$conf.int[2]
  results$p[i]       <- ft$p.value
}

results$p_fdr <- p.adjust(results$p, method = "fdr")
print(results)

# Summarize the per-region ORs
mean(results$OR)   # mean OR across regions
range(results$OR)  # range


####################
# Pooled rate
###################
cat("dogCD overall rate:", sum(dogCD_overlap) / (dogCDbp * nrow(count_matrix_dogCD)), "\n")
cat("size overall rate: ", sum(size_overlap)  / (Size_bp  * nrow(count_matrix_size)),  "\n")

# Mantel-Haenszel pooled OR across regions
library(metafor)

dogCDbp <- 6790724
Size_bp  <- 6377108

dogCD_overlap <- as.data.frame(count_matrix_dogCD)$Overlap_bp
size_overlap  <- as.data.frame(count_matrix_size)$Overlap_bp

# or simply a pooled Fisher on summed counts:
mat_pooled <- matrix(c(sum(dogCD_overlap),  dogCDbp*8 - sum(dogCD_overlap),
                       sum(size_overlap),   Size_bp*8  - sum(size_overlap)),
                     nrow = 2, byrow = TRUE)
fisher.test(mat_pooled)


# Format for mantelhaen.test: needs a 2x2xK array
arr <- array(NA, dim = c(2, 2, length(regions)),
             dimnames = list(c("dogCD","size"), c("overlap","no_overlap"), regions))

for (i in seq_along(regions)) {
  arr[,,i] <- matrix(c(dogCD_overlap[i],  dogCDbp - dogCD_overlap[i],
                       size_overlap[i],   Size_bp  - size_overlap[i]),
                     nrow = 2, byrow = TRUE)
}

mantelhaen.test(arr)

##########

library(tidyverse)
library(patchwork)

df <- tribble(
  ~region, ~rate_dogCD, ~rate_size, ~OR, ~CI_low, ~CI_high,
  "ACG", 0.06, 0.04, 1.36, 1.35, 1.37,
  "cerebellum", 0.08, 0.07, 1.14, 1.14, 1.15,
  "frontal", 0.09, 0.08, 1.17, 1.16, 1.17,
  "hypothalamus", 0.08, 0.07, 1.21, 1.21, 1.22,
  "occipital", 0.09, 0.07, 1.23, 1.22, 1.23,
  "striatum", 0.08, 0.07, 1.18, 1.17, 1.18,
  "temporal", 0.08, 0.06, 1.22, 1.22, 1.23,
  "thalamus", 0.11, 0.09, 1.22, 1.22, 1.23
) %>%
  mutate(region = factor(region, levels = rev(region)))

p1 <- df %>%
  pivot_longer(c(rate_dogCD, rate_size), names_to = "type", values_to = "rate") %>%
  ggplot(aes(x = rate, y = region, fill = type)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  scale_fill_manual(values = c("rate_dogCD" = "#1f77b4", "rate_size" = "#ff7f0e"),
                    labels = c("dogCD", "size")) +
  labs(x = "Rate", y = "Region", fill = NULL) +
  theme_minimal(base_size = 12)

p2 <- ggplot(df, aes(x = OR, y = region)) +
  geom_errorbar(
    aes(xmin = CI_low, xmax = CI_high),
    orientation = "y",
    width = 0.2
  ) +
  geom_point(size = 3) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  labs(x = "OR", y = NULL) +
  theme_minimal(base_size = 12)

p1 + p2 + plot_layout(widths = c(1, 1.1))
###
library(tidyverse)
library(patchwork)

df <- tribble(
  ~region, ~rate_dogCD, ~rate_size, ~OR, ~CI_low, ~CI_high,
  "ACG", 0.05669, 0.04237, 1.35831, 1.35143, 1.36521,
  "cerebellum", 0.08103, 0.07154, 1.14440, 1.13977, 1.14910,
  "frontal", 0.09146, 0.07941, 1.16700, 1.16245, 1.17152,
  "hypothalamus", 0.07890, 0.06602, 1.21186, 1.20677, 1.21700,
  "occipital", 0.08859, 0.07345, 1.22621, 1.22130, 1.23110,
  "striatum", 0.07790, 0.06693, 1.17784, 1.17291, 1.18278,
  "temporal", 0.07611, 0.06313, 1.22261, 1.21735, 1.22780,
  "thalamus", 0.10799, 0.09001, 1.22403, 1.21962, 1.22848
) %>%
  mutate(region = factor(region, levels = rev(region)))

p1 <- df %>%
  pivot_longer(c(rate_dogCD, rate_size), names_to = "type", values_to = "rate") %>%
  ggplot(aes(x = rate, y = region, fill = type)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  scale_fill_manual(values = c("rate_dogCD" = "#1f77b4", "rate_size" = "#ff7f0e"),
                    labels = c("dogCD", "size")) +
  labs(x = "Rate", y = "Region", fill = NULL) +
  theme_minimal(base_size = 12)

p2 <- ggplot(df, aes(x = OR, y = region)) +
  geom_errorbar(
    aes(xmin = CI_low, xmax = CI_high),
    orientation = "y",
    width = 0.2
  ) +
  geom_point(size = 3) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  labs(x = "OR", y = NULL) +
  theme_minimal(base_size = 12)

p1 + p2 + plot_layout(widths = c(1, 1.1))

###########
library(tidyverse)
library(patchwork)

df <- tribble(
  ~region, ~rate_dogCD, ~rate_size, ~OR, ~CI_low, ~CI_high,
  "ACG", 0.06, 0.04, 1.36, 1.35, 1.37,
  "cerebellum", 0.08, 0.07, 1.14, 1.14, 1.15,
  "frontal", 0.09, 0.08, 1.17, 1.16, 1.17,
  "hypothalamus", 0.08, 0.07, 1.21, 1.21, 1.22,
  "occipital", 0.09, 0.07, 1.23, 1.22, 1.23,
  "striatum", 0.08, 0.07, 1.18, 1.17, 1.18,
  "temporal", 0.08, 0.06, 1.22, 1.22, 1.23,
  "thalamus", 0.11, 0.09, 1.22, 1.22, 1.23
) %>%
  mutate(region = factor(region, levels = rev(region)))

dogCD_overlap <- c() # replace with your vector
size_overlap   <- c() # replace with your vector

dogCDbp <- 6790724
Size_bp <- 6377108

# pooled stats
dogCD_x <- sum(dogCD_overlap)
size_x  <- sum(size_overlap)
dogCD_n <- dogCDbp * length(df$region)
size_n  <- Size_bp * length(df$region)

dogCD_ci <- binom.test(dogCD_x, dogCD_n)$conf.int
size_ci  <- binom.test(size_x, size_n)$conf.int

pooled_df <- tibble(
  type = factor(c("dogCD", "size"), levels = c("size", "dogCD")),
  rate = c(dogCD_x / dogCD_n, size_x / size_n),
  low = c(dogCD_ci[1], size_ci[1]),
  high = c(dogCD_ci[2], size_ci[2]),
  n = c(dogCD_n, size_n),
  x = c(dogCD_x, size_x)
) %>%
  mutate(
    label = sprintf("%s  rate=%.3f  n=%s", type, rate, scales::comma(n))
  )

p1 <- df %>%
  pivot_longer(c(rate_dogCD, rate_size), names_to = "type", values_to = "rate") %>%
  ggplot(aes(x = rate, y = region, fill = type)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  scale_fill_manual(values = c("rate_dogCD" = "#1f77b4", "rate_size" = "#ff7f0e"),
                    labels = c("dogCD", "size")) +
  labs(x = "Rate", y = "Region", fill = NULL) +
  theme_minimal(base_size = 12)

p2 <- ggplot(df, aes(x = OR, y = region)) +
  geom_errorbar(aes(xmin = CI_low, xmax = CI_high),
                orientation = "y", width = 0.2) +
  geom_point(size = 3) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  labs(x = "OR", y = NULL) +
  theme_minimal(base_size = 12)

p3 <- ggplot(pooled_df, aes(x = rate, y = type)) +
  geom_errorbarh(aes(xmin = low, xmax = high), height = 0.18) +
  geom_point(size = 3) +
  geom_text(aes(label = label), hjust = -0.05, size = 3.2) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, max(pooled_df$high) * 1.25)) +
  labs(x = "Pooled rate", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")

p1 + p2 + p3 + plot_layout(widths = c(1.2, 1.1, 1))

######
## Forest plot: dogCD vs. body-size GWAS region overlap with brain cCREs
## Odds ratio (95% CI) per brain region, Fisher's exact test
##
## To use your own data instead of the hardcoded values below, read it
## directly (comma decimals, tab-separated, exactly as you pasted it):
##   df <- read.delim("your_table.tsv", dec = ",")

library(ggplot2)
library(dplyr)

df <- data.frame(
  region     = c("ACG", "cerebellum", "frontal", "hypothalamus",
                 "occipital", "striatum", "temporal", "thalamus"),
  rate_dogCD = c(0.06, 0.08, 0.09, 0.08, 0.09, 0.08, 0.08, 0.11),
  rate_size  = c(0.04, 0.07, 0.08, 0.07, 0.07, 0.07, 0.06, 0.09),
  OR         = c(1.36, 1.14, 1.17, 1.21, 1.23, 1.18, 1.22, 1.22),
  CI_low     = c(1.35, 1.14, 1.16, 1.21, 1.22, 1.17, 1.22, 1.22),
  CI_high    = c(1.37, 1.15, 1.17, 1.22, 1.23, 1.18, 1.23, 1.23)
)

# Order regions by OR so the strongest enrichment sits at the top
df <- df %>%
  mutate(
    region    = factor(region, levels = region[order(OR)]),
    highlight = ifelse(as.character(region) == "ACG", "ACG", "Other regions"),
    label     = sprintf("%.2f (%.2f\u2013%.2f)", OR, CI_low, CI_high)
  )

p <- ggplot(df, aes(y = region)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey60", linewidth = 0.4) +
  geom_segment(aes(x = CI_low, xend = CI_high, yend = region, color = highlight),
               linewidth = 1.0, lineend = "round") +
  geom_point(aes(x = OR, color = highlight), size = 3.2) +
  geom_text(aes(x = 1.40, label = label), hjust = 0, size = 3.1, family = "sans", color = "grey20") +
  scale_color_manual(values = c("ACG" = "#D55E00", "Other regions" = "#2C5F7C"), guide = "none") +
  scale_x_continuous(limits = c(1.00, 1.65), breaks = seq(1.0, 1.4, 0.1)) +
  labs(
    x = "Odds ratio (dogCD vs. body-size GWAS\u2013cCRE overlap)",
    y = NULL,
    title = "Enrichment of dogCD GWAS regions in brain cCREs"
  ) +
  theme_classic(base_size = 12, base_family = "sans") +
  theme(
    plot.title = element_text(size = 11, face = "plain"),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()
  )

ggsave("forest_plot_dogCD.pdf", p, width = 5.2, height = 3.6, device = cairo_pdf)
ggsave("forest_plot_dogCD.png", p, width = 5.2, height = 3.6, dpi = 300)

###
## Two-panel figure: dogCD vs. body-size GWAS region overlap with brain cCREs
## Panel A: overlap rates (dogCD vs. body size) per brain region
## Panel B: odds ratio (95% CI), Fisher's exact test, per brain region
##
## To use your own data instead of the hardcoded values below, read it
## directly (comma decimals, tab-separated, exactly as you pasted it):
##   df <- read.delim("your_table.tsv", dec = ",")

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

df <- data.frame(
  region     = c("ACG", "cerebellum", "frontal", "hypothalamus",
                 "occipital", "striatum", "temporal", "thalamus"),
  rate_dogCD = c(0.05669, 0.08103, 0.09146, 0.07890, 0.08859, 0.07790, 0.07611, 0.10799),
  rate_size  = c(0.04237, 0.07154, 0.07941, 0.06602, 0.07345, 0.06693, 0.06313, 0.09001),
  OR         = c(1.35831, 1.14440, 1.16700, 1.21186, 1.22621, 1.17784, 1.22261, 1.22403),
  CI_low     = c(1.35143, 1.13977, 1.16245, 1.20677, 1.22130, 1.17291, 1.21735, 1.21962),
  CI_high    = c(1.36521, 1.14910, 1.17152, 1.21700, 1.23110, 1.18278, 1.22780, 1.22848)
)

# Order regions by OR (ascending) so the strongest enrichment sits at the top
# in both panels, and carry that order to the rates panel for alignment
region_order <- df$region[order(df$OR)]
df <- df %>%
  mutate(
    region    = factor(region, levels = region_order),
    highlight = ifelse(as.character(region) == "ACG", "ACG", "Other regions"),
    label     = sprintf("%.3f (%.3f\u2013%.3f)", OR, CI_low, CI_high)
  )

## ---- Panel A: overlap rates ----
rates_long <- df %>%
  select(region, rate_dogCD, rate_size) %>%
  pivot_longer(cols = c(rate_dogCD, rate_size), names_to = "trait", values_to = "rate") %>%
  mutate(trait = recode(trait, rate_dogCD = "dogCD", rate_size = "Body size"))

p_rates <- ggplot(rates_long, aes(x = rate, y = region, fill = trait)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.62) +
  scale_fill_manual(values = c("dogCD" = "#D55E00", "Body size" = "grey65"), name = NULL) +
  scale_x_continuous(labels = scales::label_percent(accuracy = 1)) +
  labs(x = "Overlap rate with brain cCREs", y = NULL) +
  theme_classic(base_size = 12, base_family = "sans") +
  theme(
    legend.position = "top",
    legend.justification = "left",
    legend.margin = margin(b = -4, l = -2),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()
  )

## ---- Panel B: odds ratio forest plot ----
p_or <- ggplot(df, aes(y = region)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey60", linewidth = 0.4) +
  geom_segment(aes(x = CI_low, xend = CI_high, yend = region, color = highlight),
               linewidth = 1.0, lineend = "round") +
  geom_point(aes(x = OR, color = highlight), size = 3.2) +
  geom_text(aes(x = 1.43, label = label), hjust = 0, size = 2.95, family = "sans", color = "grey20") +
  scale_color_manual(values = c("ACG" = "#D55E00", "Other regions" = "#2C5F7C"), guide = "none") +
  scale_x_continuous(limits = c(1.00, 1.95), breaks = seq(1.0, 1.4, 0.1)) +
  labs(x = "Odds ratio (dogCD vs. body size)", y = NULL) +
  theme_classic(base_size = 12, base_family = "sans") +
  theme(
    plot.margin = margin(t = 23),  # nudge down to align with panel A under its legend
    axis.text.y = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()
  )

combined <- p_rates + p_or +
  plot_layout(widths = c(1, 1.25)) +
  plot_annotation(
    title = "dogCD GWAS region overlap with brain cCREs, by region",
    tag_levels = "A",
    theme = theme(plot.title = element_text(size = 12, family = "sans"))
  )

ggsave("forest_plot_dogCD.pdf", combined, width = 9.5, height = 4.0, device = cairo_pdf)
ggsave("forest_plot_dogCD.png", combined, width = 9.5, height = 4.0, dpi = 300)


