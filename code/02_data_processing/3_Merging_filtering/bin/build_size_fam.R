#!/usr/bin/env Rscript
# build_size_fam.R — SIZE control-phenotype fam file (run via mlma alongside F1/F2/F3)
# Reproduces only the SIZE section of data/04_gwas/getSIZE_pheno_genofiles.r; that file's second
# section (STUCK, a different Darwin's Ark item) is not reproduced here.
#
# Uses this pipeline's own already-built data1.txt directly, rather than a fresh
# famQC5+covariates join — the original instead re-reads data1 from a round trip through disk via
# a join that is otherwise identical to BUILD_DATA1's.
#
# Q121A.pheno (Darwin's Ark's own SIZE survey-derived phenotype,
# data/02_data_processing/Q121A.pheno) is exposed as the required `size_pheno_file` param.

suppressMessages(library(dplyr))

args <- commandArgs(trailingOnly = TRUE)
data1_path      <- args[1]
size_pheno_path <- args[2]

data1 <- read.csv(data1_path, header = TRUE, sep = "\t")
data1$IID <- as.character(data1$IID)

size_pheno <- read.table(size_pheno_path, header = FALSE, col.names = c("FID", "IID", "pheno"), sep = " ")
size_pheno$IID <- as.character(size_pheno$IID)

sink("build_size_fam.log", split = TRUE)

table(size_pheno$pheno, exclude = NULL)

size2 <- left_join(data1, size_pheno, by = "IID")
nrow(size2)

# Set pheno to NA where sex_numeric is NA
size2$pheno[is.na(size2$sex_numeric)] <- NA

# Recode sex_numeric: 1=1, 2=0, NA=NA
size2$sex_numeric <- ifelse(size2$sex_numeric == 1, 1, ifelse(size2$sex_numeric == 2, 0, NA))

# Recode pheno NA to -9
size2$pheno[is.na(size2$pheno)] <- -9

table(size2$sex_numeric, size2$pheno, exclude = NULL)

sink()

size_fam <- data.frame(FID = size2$IID, IID = size2$IID, F = size2$M, M = size2$M,
                        sex = size2$sex_numeric, phe = size2$pheno)

write.table(size_fam, file = "SIZE_DA_MERGED_GENCOVE_AXIOM_QC5.fam",
            row.names = FALSE, sep = "\t", col.names = FALSE, quote = FALSE)
