#!/usr/bin/env Rscript
# build_factor_fam.R — per-factor NA-filtered phenotype fam file
#
# Genericized across F1/F2/F3 via arguments. Real per-factor differences preserved:
#   - which item columns feed the NA count (F1: item7/153/154/155, F2: item93/95/150,
#     F3: item145/146/147/148/149/151/152)
#   - the NA-count exclusion threshold (F1: >=3, F2: >=2, F3: >=7) — literal per-factor values
#     from the original, not derived from item count.
#
# F1's NA-count filter uses its own `na_count_F1` column, matching F2/F3's `na_count_F{2,3}`
# usage (the original's F1 block instead references a bare, undefined `na_count`).
#
# The fam file is written to disk here (the original's equivalent write was commented out, even
# though the next step reads this exact file back in). One NA-count diagnostic PNG is produced
# per factor.

suppressMessages(library(dplyr))

args <- commandArgs(trailingOnly = TRUE)
data6_path      <- args[1]
factor_id       <- args[2]                          # "F1", "F2", "F3"
item_cols       <- strsplit(args[3], ",")[[1]]       # e.g. c("item7","item153","item154","item155")
na_threshold    <- as.numeric(args[4])

data6 <- read.csv(data6_path, sep = "\t", header = TRUE)

na_count_col <- paste0("na_count_", factor_id)
modi_col     <- paste0(factor_id, "modi")
modi2_col    <- paste0(factor_id, "modi2")
modi3_col    <- paste0(factor_id, "modi3")

sink(paste0("build_", factor_id, "_fam.log"), split = TRUE)

summary(data6[[factor_id]], exclude = NULL)
summary(data6[[paste0(factor_id, "_SE")]], exclude = NULL)

data6[[na_count_col]] <- apply(data6[, item_cols], 1, function(row) sum(is.na(row)))
table(data6[[na_count_col]], exclude = NULL)

data6[[modi_col]] <- ifelse(data6[[na_count_col]] >= na_threshold, NA, data6[[factor_id]])
table(data6[[na_count_col]], exclude = NULL)
summary(data6[[factor_id]], exclude = NULL)
summary(data6[[modi_col]], exclude = NULL)

data6[[modi2_col]] <- ifelse(is.na(data6$age), NA, data6[[modi_col]])
table(data6[[modi2_col]], exclude = NULL)
nrow(data6)

data6[[modi3_col]] <- ifelse(is.na(data6[[modi2_col]]), -9, data6[[modi2_col]])
summary(data6[[modi3_col]], exclude = NULL)
table(data6$sex, data6[[modi2_col]], exclude = NULL)

data6$M <- 0

sink()

png(paste0(factor_id, "_na_counts.png"), width = 6, height = 4, units = "in", res = 150)
plot(data6[[na_count_col]], xlab = "N", pch = 19, ylab = paste("NA count", factor_id), main = paste0("CCD", factor_id))
dev.off()

fam_out <- data.frame(FID = data6$IID, IID = data6$IID, F = data6$M, M = data6$M,
                       sex = data6$sex, phe = data6[[modi3_col]])

write.table(fam_out, file = paste0("CCD", factor_id, "_DA_MERGED_GENCOVE_AXIOM_QC5.fam"),
            row.names = FALSE, sep = "\t", col.names = FALSE, quote = FALSE)
