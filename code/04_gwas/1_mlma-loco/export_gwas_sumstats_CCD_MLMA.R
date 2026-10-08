###export the sumstats in the right format using R ####
#OBS! here GWAS sumstats from MLMA:

module load bioinfo-tools R/4.3.1
module load R_packages/4.3.1

R
setwd("/proj/snic2022-6-229/private/Darwin_Affy/NEW_DA_IMP_MERGED/2024-10-11/")

#CCDF1
Wood <- read.table("CCDF1_QC6_LOCO.loco.mlma", header=T) 
Wood$N <- 2429
Wood1 <- data.frame(SNP=Wood$SNP, A1=Wood$A1, A2=Wood$A2, freq=Wood$Freq, BETA=Wood$b, se=Wood$se, P=Wood$p, N=Wood$N, CHR=Wood$Chr, POS=Wood$bp)
write.table(Wood1, file="/proj/snic2022-6-229/private/Darwin_Affy/NEW_DA_IMP_MERGED/dog-finemapping_polyfun/CCD_250822/geno_files/CCDF1_sumstats.txt", col.names=T, row.names=F, quote=F)

#CCDF2
Wood <- read.table("CCDF2_QC6_LOCO.loco.mlma", header=T) 
Wood$N <- 2506
Wood1 <- data.frame(SNP=Wood$SNP, A1=Wood$A1, A2=Wood$A2, freq=Wood$Freq, BETA=Wood$b, se=Wood$se, P=Wood$p, N=Wood$N, CHR=Wood$Chr, POS=Wood$bp)
write.table(Wood1, file="/proj/snic2022-6-229/private/Darwin_Affy/NEW_DA_IMP_MERGED/dog-finemapping_polyfun/CCD_250822/geno_files/CCDF2_sumstats.txt", col.names=T, row.names=F, quote=F)

#CCDF3
Wood <- read.table("CCDF3_QC6_LOCO.loco.mlma", header=T) 
Wood$N <- 2428
Wood1 <- data.frame(SNP=Wood$SNP, A1=Wood$A1, A2=Wood$A2, freq=Wood$Freq, BETA=Wood$b, se=Wood$se, P=Wood$p, N=Wood$N, CHR=Wood$Chr, POS=Wood$bp)
write.table(Wood1, file="/proj/snic2022-6-229/private/Darwin_Affy/NEW_DA_IMP_MERGED/dog-finemapping_polyfun/CCD_250822/geno_files/CCDF3_sumstats.txt", col.names=T, row.names=F, quote=F)

