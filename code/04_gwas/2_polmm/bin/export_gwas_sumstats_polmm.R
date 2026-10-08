#!/usr/bin/env Rscript
# export_gwas_sumstats_polmm.R — reformat POLMM's reformatted marker output for finemapping
# Source: 04_gwas/2_polmm/export_gwas_sumstats_CCD_POLMM.R
#
# Genericized across all 14 items via arguments — the original repeats this block once per item,
# identical apart from the file paths.
#
# The A1/A2 swap (`A1 = wood$A2, A2 = wood$A1`) is intentional, not a bug: GRAB/POLMM's own A1/A2
# convention doesn't match what this sumstats format (and downstream finemapping) expects BETA to
# be measured against, so the columns are deliberately relabeled here. GRAB.Marker's output labels
# its Info field "CHR:POS:REF:ALT" and reports AltFreq/AltCounts/beta/seBeta relative to ALT, but
# reformat_polmm_output.sh's positional split labels REF as "A1" and ALT as "A2" (matching Info's
# REF-then-ALT order literally, not by meaning) — so before this swap, freq/BETA are actually
# measured against A2, not A1. The swap makes the final A1 the ALT allele, consistent with
# freq/BETA both being reported against it.
#
# N (= AltCounts/AltFreq) is computed from real per-marker values.

args <- commandArgs(trailingOnly = TRUE)
modi_path <- args[1]
out_path  <- args[2]

wood <- read.table(modi_path, header = TRUE)
wood$SNP <- paste(wood$CHR, wood$Position, sep = ":")
wood$N <- wood$AltCounts / wood$AltFreq

sumstats <- data.frame(SNP = wood$SNP, A1 = wood$A2, A2 = wood$A1, freq = wood$AltFreq,
                        BETA = wood$beta, se = wood$seBeta, P = wood$Pvalue, N = wood$N,
                        CHR = wood$Chr, POS = wood$bp)

write.table(sumstats, file = out_path, col.names = TRUE, row.names = FALSE, quote = FALSE)
