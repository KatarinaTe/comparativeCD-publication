#!/usr/bin/env Rscript
# run_polmm_item.R — per-item POLMM null model fit + full-genome marker test
# Source: 04_gwas/2_polmm/POLMMgrab_PER_ITEM.sh
#
# Genericized across all 14 items via the item_id argument. Step 1 fits the null model on the
# pruned genotype (EigenInput, ~1M SNPs); step 2 runs the full-genome marker test using that
# fitted model.

suppressMessages(library(GRAB))

args <- commandArgs(trailingOnly = TRUE)
grab_file  <- args[1]
full_bed   <- args[2]  # CCDitem{N}_DA_MERGED_GENCOVE_AXIOM_QC6.bed
pruned_bed <- args[3]  # CCDitem{N}_EigenInput.bed
item_id    <- args[4]
out_file   <- args[5]

pheno_data <- read.table(grab_file, header = TRUE, sep = "")

item_col   <- paste0("item", item_id)
factor_col <- paste0(item_col, "_factor")
pheno_data[[factor_col]] <- as.factor(pheno_data[[item_col]])

null_model_formula <- as.formula(paste0(factor_col, " ~ sex + age"))

obj_polmm <- GRAB.NullModel(formula = null_model_formula,
                             data = pheno_data,
                             subjData = pheno_data$IID,
                             method = "POLMM",
                             traitType = "ordinal",
                             GenoFile = pruned_bed,
                             control = list(showInfo = FALSE,
                                            LOCO = FALSE,
                                            tolTau = 0.2,
                                            tolBeta = 0.1))

GRAB.Marker(obj_polmm, GenoFile = full_bed, OutputFile = out_file)
