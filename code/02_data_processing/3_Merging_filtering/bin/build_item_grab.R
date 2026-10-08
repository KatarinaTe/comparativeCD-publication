#!/usr/bin/env Rscript
# build_item_grab.R — per-item combined phenotype/covariate "grab" file for POLMM
#
# `sex.x` is dplyr's automatic collision suffix: both famQC6 (from the plink fam file) and data6
# have a `sex` column, so the join produces `sex.x` (famQC6's own plink-coded sex) and `sex.y`
# (data6's sex); `sex.x` is used here.

suppressMessages(library(dplyr))

args <- commandArgs(trailingOnly = TRUE)
qc6_fam_path <- args[1]
data6_path   <- args[2]
item_id      <- args[3]   # e.g. "155"

data6 <- read.csv(data6_path, sep = "\t", header = TRUE)

item_col <- paste0("item", item_id)

famQC6 <- read.csv(qc6_fam_path, sep = " ", header = FALSE, col.names = c("FID", "IID", "F", "M", "sex", "pheno"))

sink(paste0("build_item", item_id, "_grab.log"), split = TRUE)

table(famQC6$pheno, exclude = NULL)

sink()

famQC6$IID <- as.character(famQC6$IID)

famQC6b <- left_join(famQC6, data6, by = "IID")

famQC6c <- famQC6b %>%
  mutate(sex_binary = ifelse(sex.x == 1, 1, 0))

QC6grab <- data.frame(IID = famQC6c$IID, sex = famQC6c$sex_binary, age = famQC6c$age)
QC6grab[[item_col]] <- famQC6c[[item_col]]
QC6grab <- QC6grab[, c("IID", item_col, "sex", "age")]

sink(paste0("build_item", item_id, "_grab.log"), append = TRUE, split = TRUE)
nrow(QC6grab)
sink()

write.table(QC6grab, file = paste0("grabCCDitem", item_id, "_DA_MERGED_GENCOVE_AXIOM_QC6.txt"),
            row.names = FALSE, sep = "\t", col.names = TRUE, quote = FALSE)
