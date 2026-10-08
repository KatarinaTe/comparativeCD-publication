#!/usr/bin/env Rscript
# export_gwas_sumstats_mlma.R — reformat gcta64 --mlma-loco output for finemapping (PolyFun)
# Source: 04_gwas/1_mlma-loco/export_gwas_sumstats_CCD_MLMA.R
#
# Genericized across F1/F2/F3 via arguments — the original repeats this block 3 times,
# identical apart from the file paths and a hardcoded N (2429/2506/2428).
#
# N is computed here from the real .phen file's row count (the actual number of dogs the GWAS
# was run on) rather than a hardcoded historical value.

args <- commandArgs(trailingOnly = TRUE)
mlma_path   <- args[1]
phen_path   <- args[2]
out_path    <- args[3]

wood <- read.table(mlma_path, header = TRUE)
n_dogs <- nrow(read.table(phen_path, header = FALSE))

sumstats <- data.frame(SNP = wood$SNP, A1 = wood$A1, A2 = wood$A2, freq = wood$Freq,
                        BETA = wood$b, se = wood$se, P = wood$p, N = n_dogs,
                        CHR = wood$Chr, POS = wood$bp)

write.table(sumstats, file = out_path, col.names = TRUE, row.names = FALSE, quote = FALSE)
