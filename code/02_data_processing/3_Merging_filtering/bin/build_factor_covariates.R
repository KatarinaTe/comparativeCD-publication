#!/usr/bin/env Rscript
# build_factor_covariates.R — per-factor age qcovar + sex covar files
#
# `data1` (BUILD_DATA1's output) is joined in here, in place of the undefined `data4` the
# original's `left_join(famQC6, data4, by="IID")` references — data1 is the object in scope with
# both an `age` and a `sex_numeric` column, joining cleanly 1:1 on IID.
#
# `data1$IID` is explicitly re-cast to character before the join: `data1` is read back here with
# `read.table()` from a plain-text round trip, which does not preserve build_data1.R's own
# character cast on write, and an all-numeric IID column would otherwise risk a `left_join`
# "incompatible types" error under a modern dplyr.
#
# The qcovar/covar files are written to disk here; the original's equivalent `write.table()`
# calls for both files are commented out for all three factors.

suppressMessages(library(dplyr))

args <- commandArgs(trailingOnly = TRUE)
qc6_fam_path <- args[1]
data1_path   <- args[2]
prefix       <- args[3]   # e.g. "CCDF1_DA_MERGED_GENCOVE_AXIOM_QC6"

data1 <- read.table(data1_path, header = TRUE, sep = "\t")
data1$IID <- as.character(data1$IID)

famQC6 <- read.csv(qc6_fam_path, sep = " ", header = FALSE, col.names = c("FID", "IID", "F", "M", "sex", "pheno"))
famQC6$IID <- as.character(famQC6$IID)

famQC6age <- left_join(famQC6, data1, by = "IID")

sink(paste0("build_", prefix, "_covariates.log"), split = TRUE)

nrow(famQC6)
nrow(famQC6age)
table(famQC6age$sex_numeric, exclude = NULL)

sink()

famQC6age.qcovar <- data.frame(FID = famQC6age$IID, IID = famQC6age$IID, age = famQC6age$age)

write.table(famQC6age.qcovar, file = paste0("age", prefix, ".qcovar"),
            row.names = FALSE, sep = "\t", col.names = FALSE, quote = FALSE)

famQC6sex <- famQC6age %>%
  mutate(sex_binary = ifelse(sex_numeric == 1, 1, 0))

sink(paste0("build_", prefix, "_covariates.log"), append = TRUE, split = TRUE)
table(famQC6sex$sex_binary, exclude = NULL)
sink()

famQC6sex.covar <- data.frame(FID = famQC6sex$IID, IID = famQC6sex$IID, sex = famQC6sex$sex_binary)

write.table(famQC6sex.covar, file = paste0("sex", prefix, ".covar"),
            row.names = FALSE, sep = "\t", col.names = FALSE, quote = FALSE)
