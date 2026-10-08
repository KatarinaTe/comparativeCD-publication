#!/usr/bin/env Rscript
# build_data1.R — join the final QC5 fam with covariates -> data1
# Source: 02_data_processing/3_Merging_filtering/01_gencove_axiom_merging_filtering.sh (tail)
#
# The original writes this with `col.names=F` (commented out in the source), but the very next
# script (03_create_files_per_factor_mlma.R) reads data1.txt back with `header=T` — writing
# without a header would make that read misalign every column. Since data1.txt isn't a frozen
# file anywhere in this repo (unlike data6/response_df_dogCDitems.txt, this one is genuinely
# produced fresh here), written with a header to match its own declared consumer.

suppressMessages(library(dplyr))

args <- commandArgs(trailingOnly = TRUE)
fam_qc5_path   <- args[1]
covariates_path <- args[2]

famQC5 <- read.csv(fam_qc5_path, sep = " ", header = FALSE, col.names = c("FID", "IID", "F", "M", "sex", "phe"))
famQC5$IID <- as.character(famQC5$IID)

updated_covar <- read.csv(covariates_path, sep = "\t", header = TRUE)
updated_covar$IID <- as.character(updated_covar$dog)

data1 <- left_join(famQC5, updated_covar, by = "IID")
write.table(data1, file = "data1.txt", row.names = FALSE, sep = "\t", col.names = TRUE, quote = FALSE)
