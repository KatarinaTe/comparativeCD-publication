#!/usr/bin/env Rscript
# build_clump_regions_100kb.R — one BED region per significant GWAS clump, spanning that
# clump's full extent (lead SNP + every secondary SNP plink assigned to it) padded by 100kb on
# each side, pooled across every phenotype passed in (e.g. CCDF1+CCDF2+CCDF3 for "dogCD", or
# SIZE alone) and merged where overlapping.
#
# plink's .clumped file only contains clumps whose lead SNP passed --clump-p1 -- but that's a
# looser "suggestive" threshold (1e-6, CLUMP_FACTOR/CLUMP_ITEM's own --clump-p1) used to catch
# every candidate clump for manual review, not the genome-wide significance threshold itself.
# This script filters down to genome-wide significant clumps only (lead SNP P < GWS_P_THRESHOLD,
# 4e-7 per the manuscript's own stated threshold) -- confirmed directly with the manuscript author
# (2026-09-23) that this is the automated equivalent of how the original GWAS-region reference
# files were built (broad clumping first, manual selection of genome-wide significant clumps
# after). It gives the lead SNP's own BP directly, but only lists secondary SNPs by ID (its SP2
# column) — their own BP positions are looked up from the matching --mlma-loco summary-stats file.
#
# Usage: build_clump_regions_100kb.R <out.bed> <clumped1> <mlma1> [<clumped2> <mlma2> ...]
# One <clumped>/<mlma> pair per phenotype being pooled into this output.
#
# The second file of each pair can be either an mlma-loco summary-stats file (native SNP + bp
# columns) or a POLMM item's own reformatted marker output (modi_simuMarkerOutput_POLMM_item*,
# see 04_gwas/2_polmm/bin/reformat_polmm_output.sh) -- those have CHR/Position/bp but no SNP
# column, since GRAB.Marker never computes one. SNP is constructed below when missing, as
# CHR:Position, matching the upstream `bcftools annotate --set-id '%CHROM:%POS'` convention that
# plink's own SP2 secondary-SNP IDs are drawn from -- confirmed this is what CLUMP_ITEM actually
# clumped against (2026-09-23, added when 06b-clump-regions-100kb was extended to include POLMM).

args <- commandArgs(trailingOnly = TRUE)
out_path  <- args[1]
pair_args <- args[-1]
stopifnot("clumped/mlma files must be given in pairs" = length(pair_args) %% 2 == 0)

flank <- 100000
gws_p_threshold <- 4e-7

region_rows <- list()

for (i in seq(1, length(pair_args), by = 2)) {
  clumped_path <- pair_args[i]
  mlma_path    <- pair_args[i + 1]

  clumped <- read.table(clumped_path, header = TRUE, stringsAsFactors = FALSE)
  clumped <- clumped[clumped$P < gws_p_threshold, ]
  if (nrow(clumped) == 0) next

  mlma <- read.table(mlma_path, header = TRUE, stringsAsFactors = FALSE)
  if (is.null(mlma$SNP)) mlma$SNP <- paste(mlma$CHR, mlma$Position, sep = ":")
  bp_by_snp <- setNames(mlma$bp, mlma$SNP)

  for (j in seq_len(nrow(clumped))) {
    sp2 <- clumped$SP2[j]
    secondary_ids <- character(0)
    if (!is.na(sp2) && sp2 != "NONE") {
      secondary_ids <- trimws(gsub("\\([^)]*\\)", "", strsplit(sp2, ",")[[1]]))
    }
    secondary_bp <- unname(bp_by_snp[secondary_ids])
    secondary_bp <- secondary_bp[!is.na(secondary_bp)]

    clump_bp <- c(clumped$BP[j], secondary_bp)

    region_rows[[length(region_rows) + 1]] <- data.frame(
      chr   = paste0("chr", clumped$CHR[j]),
      start = max(0, min(clump_bp) - flank),
      end   = max(clump_bp) + flank
    )
  }
}

stopifnot("no significant clumps found in any input file" = length(region_rows) > 0)

regions <- do.call(rbind, region_rows)
regions <- regions[order(regions$chr, regions$start), ]

write.table(regions, out_path, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
