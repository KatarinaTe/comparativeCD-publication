##export the sumstats in the right format using R ####

#OBS! here GWAS sumstats from POLMM:
module load bioinfo-tools R/4.3.1
module load R_packages/4.3.1

R
setwd("/proj/snic2022-6-229/private/Darwin_Affy/NEW_DA_IMP_MERGED/2024-10-11/")


#item145
Wood <- read.table("modi_simuMarkerOutput_POLMM_item145_FULLGRM_fullGeno.txt", header=T)
Wood$SNP <- paste(Wood$CHR, Wood$Position, sep = ":")
Wood$N <- Wood$AltCounts/Wood$AltFreq
#OBS! BETA must be for the A1 allele!!!!!!!!!!!!!!!!!!!!!!!!! #changing the A2 to become A1:
Wood1 <- data.frame(SNP=Wood$SNP, A1=Wood$A2, A2=Wood$A1, freq=Wood$AltFreq, BETA=Wood$beta, se=Wood$seBeta, P=Wood$Pvalue, N=Wood$N, CHR=Wood$Chr, POS=Wood$bp)
write.table(Wood1, file="/proj/snic2022-6-229/private/Darwin_Affy/NEW_DA_IMP_MERGED/dog-finemapping_polyfun/CCD_250822/ITEM145_sumstats.txt", col.names=T, row.names=F, quote=F)

#item146
Wood <- read.table("modi_simuMarkerOutput_POLMM_item146_FULLGRM_fullGeno.txt", header=T)
Wood$SNP <- paste(Wood$CHR, Wood$Position, sep = ":")
Wood$N <- Wood$AltCounts/Wood$AltFreq
#OBS! BETA must be for the A1 allele!!!!!!!!!!!!!!!!!!!!!!!!! #changing the A2 to become A1:
Wood1 <- data.frame(SNP=Wood$SNP, A1=Wood$A2, A2=Wood$A1, freq=Wood$AltFreq, BETA=Wood$beta, se=Wood$seBeta, P=Wood$Pvalue, N=Wood$N, CHR=Wood$Chr, POS=Wood$bp)
write.table(Wood1, file="/proj/snic2022-6-229/private/Darwin_Affy/NEW_DA_IMP_MERGED/dog-finemapping_polyfun/CCD_250822/geno_files/ITEM146_sumstats.txt", col.names=T, row.names=F, quote=F)

#item147
Wood <- read.table("modi_simuMarkerOutput_POLMM_item147_FULLGRM_fullGeno.txt", header=T)
Wood$SNP <- paste(Wood$CHR, Wood$Position, sep = ":")
Wood$N <- Wood$AltCounts/Wood$AltFreq
#OBS! BETA must be for the A1 allele!!!!!!!!!!!!!!!!!!!!!!!!! #changing the A2 to become A1:
Wood1 <- data.frame(SNP=Wood$SNP, A1=Wood$A2, A2=Wood$A1, freq=Wood$AltFreq, BETA=Wood$beta, se=Wood$seBeta, P=Wood$Pvalue, N=Wood$N, CHR=Wood$Chr, POS=Wood$bp)
write.table(Wood1, file="/proj/snic2022-6-229/private/Darwin_Affy/NEW_DA_IMP_MERGED/dog-finemapping_polyfun/CCD_250822/geno_files/ITEM147_sumstats.txt", col.names=T, row.names=F, quote=F)

#item148
Wood <- read.table("modi_simuMarkerOutput_POLMM_item148_FULLGRM_fullGeno.txt", header=T)
Wood$SNP <- paste(Wood$CHR, Wood$Position, sep = ":")
Wood$N <- Wood$AltCounts/Wood$AltFreq
#OBS! BETA must be for the A1 allele!!!!!!!!!!!!!!!!!!!!!!!!! #changing the A2 to become A1:
Wood1 <- data.frame(SNP=Wood$SNP, A1=Wood$A2, A2=Wood$A1, freq=Wood$AltFreq, BETA=Wood$beta, se=Wood$seBeta, P=Wood$Pvalue, N=Wood$N, CHR=Wood$Chr, POS=Wood$bp)
write.table(Wood1, file="/proj/snic2022-6-229/private/Darwin_Affy/NEW_DA_IMP_MERGED/dog-finemapping_polyfun/CCD_250822/geno_files/ITEM148_sumstats.txt", col.names=T, row.names=F, quote=F)

#item149
Wood <- read.table("modi_simuMarkerOutput_POLMM_item149_FULLGRM_fullGeno.txt", header=T)
Wood$SNP <- paste(Wood$CHR, Wood$Position, sep = ":")
Wood$N <- Wood$AltCounts/Wood$AltFreq
#OBS! BETA must be for the A1 allele!!!!!!!!!!!!!!!!!!!!!!!!! #changing the A2 to become A1:
Wood1 <- data.frame(SNP=Wood$SNP, A1=Wood$A2, A2=Wood$A1, freq=Wood$AltFreq, BETA=Wood$beta, se=Wood$seBeta, P=Wood$Pvalue, N=Wood$N, CHR=Wood$Chr, POS=Wood$bp)
write.table(Wood1, file="/proj/snic2022-6-229/private/Darwin_Affy/NEW_DA_IMP_MERGED/dog-finemapping_polyfun/CCD_250822/geno_files/ITEM149_sumstats.txt", col.names=T, row.names=F, quote=F)

#item150
Wood <- read.table("modi_simuMarkerOutput_POLMM_item150_FULLGRM_fullGeno.txt", header=T)
Wood$SNP <- paste(Wood$CHR, Wood$Position, sep = ":")
Wood$N <- Wood$AltCounts/Wood$AltFreq
#OBS! BETA must be for the A1 allele!!!!!!!!!!!!!!!!!!!!!!!!! #changing the A2 to become A1:
Wood1 <- data.frame(SNP=Wood$SNP, A1=Wood$A2, A2=Wood$A1, freq=Wood$AltFreq, BETA=Wood$beta, se=Wood$seBeta, P=Wood$Pvalue, N=Wood$N, CHR=Wood$Chr, POS=Wood$bp)
write.table(Wood1, file="/proj/snic2022-6-229/private/Darwin_Affy/NEW_DA_IMP_MERGED/dog-finemapping_polyfun/CCD_250822/geno_files/ITEM150_sumstats.txt", col.names=T, row.names=F, quote=F)

#item151
Wood <- read.table("modi_simuMarkerOutput_POLMM_item151_FULLGRM_fullGeno.txt", header=T)
Wood$SNP <- paste(Wood$CHR, Wood$Position, sep = ":")
Wood$N <- Wood$AltCounts/Wood$AltFreq
#OBS! BETA must be for the A1 allele!!!!!!!!!!!!!!!!!!!!!!!!! #changing the A2 to become A1:
Wood1 <- data.frame(SNP=Wood$SNP, A1=Wood$A2, A2=Wood$A1, freq=Wood$AltFreq, BETA=Wood$beta, se=Wood$seBeta, P=Wood$Pvalue, N=Wood$N, CHR=Wood$Chr, POS=Wood$bp)
write.table(Wood1, file="/proj/snic2022-6-229/private/Darwin_Affy/NEW_DA_IMP_MERGED/dog-finemapping_polyfun/CCD_250822/geno_files/ITEM151_sumstats.txt", col.names=T, row.names=F, quote=F)

#item152
Wood <- read.table("modi_simuMarkerOutput_POLMM_item152_FULLGRM_fullGeno.txt", header=T)
Wood$SNP <- paste(Wood$CHR, Wood$Position, sep = ":")
Wood$N <- Wood$AltCounts/Wood$AltFreq
#OBS! BETA must be for the A1 allele!!!!!!!!!!!!!!!!!!!!!!!!! #changing the A2 to become A1:
Wood1 <- data.frame(SNP=Wood$SNP, A1=Wood$A2, A2=Wood$A1, freq=Wood$AltFreq, BETA=Wood$beta, se=Wood$seBeta, P=Wood$Pvalue, N=Wood$N, CHR=Wood$Chr, POS=Wood$bp)
write.table(Wood1, file="/proj/snic2022-6-229/private/Darwin_Affy/NEW_DA_IMP_MERGED/dog-finemapping_polyfun/CCD_250822/geno_files/ITEM152_sumstats.txt", col.names=T, row.names=F, quote=F)

#item153
Wood <- read.table("modi_simuMarkerOutput_POLMM_item153_FULLGRM_fullGeno.txt", header=T)
Wood$SNP <- paste(Wood$CHR, Wood$Position, sep = ":")
Wood$N <- Wood$AltCounts/Wood$AltFreq
#OBS! BETA must be for the A1 allele!!!!!!!!!!!!!!!!!!!!!!!!! #changing the A2 to become A1:
Wood1 <- data.frame(SNP=Wood$SNP, A1=Wood$A2, A2=Wood$A1, freq=Wood$AltFreq, BETA=Wood$beta, se=Wood$seBeta, P=Wood$Pvalue, N=Wood$N, CHR=Wood$Chr, POS=Wood$bp)
write.table(Wood1, file="/proj/snic2022-6-229/private/Darwin_Affy/NEW_DA_IMP_MERGED/dog-finemapping_polyfun/CCD_250822/geno_files/ITEM153_sumstats.txt", col.names=T, row.names=F, quote=F)

#item154
Wood <- read.table("modi_simuMarkerOutput_POLMM_item154_FULLGRM_fullGeno.txt", header=T)
Wood$SNP <- paste(Wood$CHR, Wood$Position, sep = ":")
Wood$N <- Wood$AltCounts/Wood$AltFreq
#OBS! BETA must be for the A1 allele!!!!!!!!!!!!!!!!!!!!!!!!! #changing the A2 to become A1:
Wood1 <- data.frame(SNP=Wood$SNP, A1=Wood$A2, A2=Wood$A1, freq=Wood$AltFreq, BETA=Wood$beta, se=Wood$seBeta, P=Wood$Pvalue, N=Wood$N, CHR=Wood$Chr, POS=Wood$bp)
write.table(Wood1, file="/proj/snic2022-6-229/private/Darwin_Affy/NEW_DA_IMP_MERGED/dog-finemapping_polyfun/CCD_250822/geno_files/ITEM154_sumstats.txt", col.names=T, row.names=F, quote=F)

#item155
Wood <- read.table("modi_simuMarkerOutput_POLMM_item155_FULLGRM_fullGeno.txt", header=T)
Wood$SNP <- paste(Wood$CHR, Wood$Position, sep = ":")
Wood$N <- Wood$AltCounts/Wood$AltFreq
#OBS! BETA must be for the A1 allele!!!!!!!!!!!!!!!!!!!!!!!!! #changing the A2 to become A1:
Wood1 <- data.frame(SNP=Wood$SNP, A1=Wood$A2, A2=Wood$A1, freq=Wood$AltFreq, BETA=Wood$beta, se=Wood$seBeta, P=Wood$Pvalue, N=Wood$N, CHR=Wood$Chr, POS=Wood$bp)
write.table(Wood1, file="/proj/snic2022-6-229/private/Darwin_Affy/NEW_DA_IMP_MERGED/dog-finemapping_polyfun/CCD_250822/geno_files/ITEM155_sumstats.txt", col.names=T, row.names=F, quote=F)

#item7
Wood <- read.table("modi_simuMarkerOutput_POLMM_item7_FULLGRM_fullGeno.txt", header=T)
Wood$SNP <- paste(Wood$CHR, Wood$Position, sep = ":")
Wood$N <- Wood$AltCounts/Wood$AltFreq
#OBS! BETA must be for the A1 allele!!!!!!!!!!!!!!!!!!!!!!!!! #changing the A2 to become A1:
Wood1 <- data.frame(SNP=Wood$SNP, A1=Wood$A2, A2=Wood$A1, freq=Wood$AltFreq, BETA=Wood$beta, se=Wood$seBeta, P=Wood$Pvalue, N=Wood$N, CHR=Wood$Chr, POS=Wood$bp)
write.table(Wood1, file="/proj/snic2022-6-229/private/Darwin_Affy/NEW_DA_IMP_MERGED/dog-finemapping_polyfun/CCD_250822/geno_files/ITEM7_sumstats.txt", col.names=T, row.names=F, quote=F)

#item93
Wood <- read.table("modi_simuMarkerOutput_POLMM_item93_FULLGRM_fullGeno.txt", header=T)
Wood$SNP <- paste(Wood$CHR, Wood$Position, sep = ":")
Wood$N <- Wood$AltCounts/Wood$AltFreq
#OBS! BETA must be for the A1 allele!!!!!!!!!!!!!!!!!!!!!!!!! #changing the A2 to become A1:
Wood1 <- data.frame(SNP=Wood$SNP, A1=Wood$A2, A2=Wood$A1, freq=Wood$AltFreq, BETA=Wood$beta, se=Wood$seBeta, P=Wood$Pvalue, N=Wood$N, CHR=Wood$Chr, POS=Wood$bp)
write.table(Wood1, file="/proj/snic2022-6-229/private/Darwin_Affy/NEW_DA_IMP_MERGED/dog-finemapping_polyfun/CCD_250822/geno_files/ITEM93_sumstats.txt", col.names=T, row.names=F, quote=F)

#item95
Wood <- read.table("modi_simuMarkerOutput_POLMM_item95_FULLGRM_fullGeno.txt", header=T)
Wood$SNP <- paste(Wood$CHR, Wood$Position, sep = ":")
Wood$N <- Wood$AltCounts/Wood$AltFreq
#OBS! BETA must be for the A1 allele!!!!!!!!!!!!!!!!!!!!!!!!! #changing the A2 to become A1:
Wood1 <- data.frame(SNP=Wood$SNP, A1=Wood$A2, A2=Wood$A1, freq=Wood$AltFreq, BETA=Wood$beta, se=Wood$seBeta, P=Wood$Pvalue, N=Wood$N, CHR=Wood$Chr, POS=Wood$bp)
write.table(Wood1, file="/proj/snic2022-6-229/private/Darwin_Affy/NEW_DA_IMP_MERGED/dog-finemapping_polyfun/CCD_250822/geno_files/ITEM95_sumstats.txt", col.names=T, row.names=F, quote=F)