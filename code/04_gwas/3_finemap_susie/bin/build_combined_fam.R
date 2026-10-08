#!/usr/bin/env Rscript
# build_combined_fam.R — builds the combined-sample-inclusion fam used to derive the one shared
# genofile that all phenotypes' finemapping runs are performed against.
# Source: 04_gwas/3_finemap_susie/1_finemap_create_genofiles.sh, lines 6-18
#
# Reconstructed using 3_Merging_filtering's own plain, whole-cohort published QC5 fam
# (DA_MERGED_GENCOVE_AXIOM_QC5.fam) rather than one of that stage's several phenotype-specific QC5
# fam variants (SIZE_..._QC5.fam, CCD${factor}_..._QC5.fam, CCDitem${id}_..._QC5.fam — same QC5
# sample set, each with a different phenotype column substituted in). This is the only one of them
# that makes sense here: this script builds ONE shared genofile used by every phenotype's
# finemapping, so a phenotype-specific starting fam would be an arbitrary pick among equally-
# eligible candidates — and whichever variant was used, only its FID/IID/F/M/sex columns would
# matter anyway, since the phenotype column gets overwritten immediately below (all QC5-stage fam
# variants describe the same underlying sample set, just with different phenotype values).
#
# Any dog present in the SIZE, STUCK, or ALLFAM phenotype fams (three fixed, already-deposited data
# files, not built by any pipeline) is kept (phenotype column set to 1); everyone else gets plink's
# missing-phenotype code (-9). This is purely a sample-inclusion mask for the shared genofile — it
# is not itself a GWAS phenotype, so STUCK's appearance here does not reintroduce STUCK as an
# analysed trait (STUCK was excluded from GWAS entirely per 3_Merging_filtering's own scope).

args <- commandArgs(trailingOnly = TRUE)
qc5_fam_path    <- args[1]
size_fam_path   <- args[2]
stuck_fam_path  <- args[3]
allfam_fam_path <- args[4]
out_fam_path    <- args[5]

fam_cols <- c("FID", "IID", "F", "M", "sex", "phe")

fam_qc5 <- read.table(qc5_fam_path, header = FALSE, col.names = fam_cols)
size    <- read.table(size_fam_path, header = FALSE, col.names = fam_cols)
stuck   <- read.table(stuck_fam_path, header = FALSE, col.names = fam_cols)
ccd     <- read.table(allfam_fam_path, header = FALSE, col.names = fam_cols)

all_IIDs <- unique(c(size$IID, ccd$IID, stuck$IID))
fam_qc5$in_any <- ifelse(fam_qc5$IID %in% all_IIDs, 1, -9)

fam_qc5$phe <- fam_qc5$in_any
fam_qc5$in_any <- NULL

write.table(fam_qc5, out_fam_path, quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)
