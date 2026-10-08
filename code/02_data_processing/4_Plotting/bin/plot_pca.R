#!/usr/bin/env Rscript
# plot_pca.R — PCA scatter plots, colored/shaped by breed
#
# `dog_breeds_csv` (data/02_data_processing/DarwinsArk_20220715_dogs_genotyped_breed_sex.csv) is a required param;
# `dog`, `purebred`, `breed1_inputted` and `sex` columns are used (`sex` appears as `sex.y` after the join).
#
# The axis labels' variance percentages ("PCA1 (Variance = 37.3%)") are computed from the real
# eigenval file rather than hardcoded.
#
# Only the colored/shaped-by-breed plot is produced. The original's first plot (uncolored, pch
# only) is not reproduced: its legend listed breed colors that never appeared on the actual
# uncolored plot.

suppressMessages(library(dplyr))

args <- commandArgs(trailingOnly = TRUE)
eigenvec_path    <- args[1]
eigenval_path    <- args[2]
data6_path       <- args[3]
dog_breeds_path  <- args[4]

pca <- read.csv(eigenvec_path, sep = " ", header = FALSE,
                 col.names = c("FID", "IID", paste0("PC", 1:20)))
pca$IID <- as.character(pca$IID)

data6 <- read.csv(data6_path, sep = "\t", header = TRUE)
data6$IID <- as.character(data6$IID)

fam_pc <- full_join(pca, data6, by = "IID")

eigenval <- read.csv(eigenval_path, sep = " ", header = FALSE, col.names = c("eigenvals"))

dog <- read.csv(dog_breeds_path, sep = ",")
dog$IID <- as.character(dog$dog)

fam_breed <- left_join(fam_pc, dog, by = "IID")

not_purebred <- subset(fam_breed, fam_breed$purebred == "no")
purebred     <- subset(fam_breed, fam_breed$purebred == "yes")

sink("pca_diagnostics.log", split = TRUE)

head(pca)
nrow(pca)
nrow(fam_pc)
head(fam_pc)
head(eigenval)
head(dog)
nrow(dog)
length(unique(dog$dog))
head(fam_breed)
table(fam_breed$breed1_inputted, exclude = NULL)
unique(fam_breed$breed1_inputted, exclude = NULL)
table(fam_breed$purebred, fam_breed$sex.y, exclude = NULL)
table(fam_breed$purebred, fam_breed$breed1_inputted, exclude = NULL)
nrow(fam_breed)
table(not_purebred$breed1_inputted, exclude = NULL)
nrow(not_purebred)
unique(purebred$breed1_inputted, exclude = NULL)
nrow(purebred)
table(purebred$sex.y, exclude = NULL)
table(fam_breed$breed1_inputted, fam_breed$purebred, exclude = NULL)

sink()

# Variance explained by PC1/PC2 — computed from the real eigenval file (see header comment).
variance_pct <- 100 * eigenval$eigenvals / sum(eigenval$eigenvals)
pc1_label <- sprintf("PCA1 (Variance = %.1f%%)", variance_pct[1])
pc2_label <- sprintf("PCA2 (Variance = %.1f%%)", variance_pct[2])

# ── Colored + shaped by breed/purebred status ───────────────────────────────────────────────────

col.vector <- rep(1, nrow(fam_breed))
col.vector[fam_breed$purebred == "yes"] <- "black"
col.vector[is.na(fam_breed$purebred)] <- "lightgrey"
col.vector[fam_breed$purebred == "no"] <- "lightgrey"
col.vector[fam_breed$breed1_inputted == "Labrador Retriever" & fam_breed$purebred == "yes"] <- "red"
col.vector[fam_breed$breed1_inputted == "Golden Retriever" & fam_breed$purebred == "yes"] <- "blue"
col.vector[fam_breed$breed1_inputted == "German Shepherd Dog" & fam_breed$purebred == "yes"] <- "purple"
col.vector[fam_breed$breed1_inputted == "Border Collie" & fam_breed$purebred == "yes"] <- "green"
col.vector[fam_breed$breed1_inputted == "Australian Shepherd" & fam_breed$purebred == "yes"] <- "yellow"
col.vector[fam_breed$breed1_inputted == "American Pit Bull Terrier" & fam_breed$purebred == "yes"] <- "orange"
col.vector[fam_breed$breed1_inputted == "Boxer" & fam_breed$purebred == "yes"] <- "lightblue3"
col.vector[fam_breed$breed1_inputted == "Siberian Husky" & fam_breed$purebred == "yes"] <- "darkred"
col.vector[fam_breed$breed1_inputted == "Chow Chow" & fam_breed$purebred == "yes"] <- "red3"
col.vector[fam_breed$breed1_inputted == "Shiba Inu" & fam_breed$purebred == "yes"] <- "green4"
col.vector[fam_breed$breed1_inputted == "Wolf"] <- "black"
col.vector[fam_breed$breed1_inputted == "New Guinea Singing Dog"] <- "yellow3"

pch.vector <- rep(1, nrow(fam_breed))
pch.vector[fam_breed$purebred == "no"] <- 2
pch.vector[fam_breed$purebred == "yes"] <- 17
pch.vector[is.na(fam_breed$purebred)] <- 21
pch.vector[fam_breed$breed1_inputted == "Wolf"] <- 8

png("pca_by_breed.png", width = 7, height = 7, units = "in", res = 150)
plot(fam_breed$PC1, fam_breed$PC2, col = col.vector, pch = pch.vector,
     xlab = pc1_label, ylab = pc2_label)
legend("topright", c("other 'purebreed' ", "Labrador Retriever", "Golden Retriever", "German Shepherd Dog", "Border Collie",
                      "Australian Shepherd", "American Pit Bull Terrier", "Boxer", "Siberian Husky", "Chow Chow", "Shiba Inu",
                      "New Guinea Singing Dog"),
       col = c("black", "red", "blue", "purple", "green", "yellow", "orange", "lightblue3", "darkred", "red3", "green4", "yellow3"),
       pch = 17, cex = 0.8)
legend("bottomright", c("purebred", "not purebred", "purebred=NA", "wolf"), pch = c(17, 2, 21, 8), cex = 0.8)
dev.off()
