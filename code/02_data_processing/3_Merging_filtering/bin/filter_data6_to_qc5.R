#!/usr/bin/env Rscript
# filter_data6_to_qc5.R — restrict + reorder data6 to exactly QC5's dogID set, in QC5's own
# row order
#
# BUILD_FACTOR_FAM/BUILD_ITEM_FAM build their per-factor/per-item fam files directly from data6,
# one row per data6 row, which FILTER_FACTOR_QC6/FILTER_ITEM_QC6 then substitute into QC5's
# .bed/.bim via plink's strictly *positional* --fam substitution. data6 is a frozen checkpoint,
# untouched by this pipeline's own genotype-side sample exclusions (LowPass-12, the dog-3094
# duplicate pair) — it still has all 3328 rows of the original QC5 population, while QC5 itself
# is now 3316 rows after those exclusions. Left unreconciled, every downstream --fam substitution
# would fail: the .bed's byte size is fixed by QC5's 3316 individuals, and a 3328-row --fam file
# can't be loaded against it.
#
# Restricts data6 to exactly QC5's own fam's dogID list, in QC5's own row order (via match(), not
# a plain ID filter — order matters as much as membership for a positional substitution). Errors
# loudly if any QC5 dog is missing from data6, rather than filtering leniently.

args <- commandArgs(trailingOnly = TRUE)
data6_path   <- args[1]
qc5_fam_path <- args[2]

data6 <- read.csv(data6_path, sep = "\t", header = TRUE)
data6$IID <- as.character(data6$IID)

qc5_fam <- read.csv(qc5_fam_path, sep = " ", header = FALSE,
                     col.names = c("FID", "IID", "F", "M", "sex", "pheno"))
qc5_fam$IID <- as.character(qc5_fam$IID)

idx <- match(qc5_fam$IID, data6$IID)
if (any(is.na(idx))) {
  stop("QC5 dogIDs missing from data6: ", paste(qc5_fam$IID[is.na(idx)], collapse = ", "))
}

data6_filtered <- data6[idx, ]

write.table(data6_filtered, file = "data6_filtered.txt",
            row.names = FALSE, sep = "\t", quote = FALSE)
