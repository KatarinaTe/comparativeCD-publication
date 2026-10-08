#!/usr/bin/env Rscript
# build_item_fam.R — per-item age-NA-filtered phenotype fam file, genericized across all 14
# items via the `item_id` argument.
#
# Unlike the per-factor fam files, there is no NA-count/multi-item threshold step here: each item
# is filtered only on age-NA, since these are single-item ordinal scores, not combined factor
# scores.
#
# Note: this script's NA counts on the committed data6_2024-10-14.txt checkpoint do not match the
# original script's own inline comment counts for item155 (728 vs. 629 raw NA, 1006 vs. 913 after
# the age-NA filter) — see archive/conversion_notes/02_data_processing.md for the full
# investigation.

suppressMessages(library(dplyr))

args <- commandArgs(trailingOnly = TRUE)
data6_path <- args[1]
item_id    <- args[2]   # e.g. "155"

data6 <- read.csv(data6_path, sep = "\t", header = TRUE)

item_col  <- paste0("item", item_id)
modi_col  <- paste0(item_col, "_modi")
modi2_col <- paste0(item_col, "_modi2")

sink(paste0("build_item", item_id, "_fam.log"), split = TRUE)

table(data6[[item_col]], exclude = NULL)

data6[[modi_col]] <- ifelse(is.na(data6$age), NA, data6[[item_col]])
table(data6[[modi_col]], exclude = NULL)
nrow(data6)
table(data6[[modi_col]], data6$sex, exclude = NULL)
# no need to exclude dogs without sex info, they don't have the phenotype anyway

data6[[modi2_col]] <- ifelse(is.na(data6[[modi_col]]), -9, data6[[modi_col]])
table(data6[[modi2_col]])

data6$M <- 0

sink()

fam_out <- data.frame(FID = data6$IID, IID = data6$IID, F = data6$M, M = data6$M,
                       sex = data6$sex, phe = data6[[modi2_col]])

write.table(fam_out, file = paste0("CCDitem", item_id, "_DA_MERGED_GENCOVE_AXIOM_QC5.fam"),
            row.names = FALSE, sep = "\t", col.names = FALSE, quote = FALSE)
