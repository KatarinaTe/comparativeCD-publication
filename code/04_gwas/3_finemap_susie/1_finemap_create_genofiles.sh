###################################################################
###these scripts use these geno inputs:
#CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.$chr
#created like these:

size <- read.table("SIZE_DA_MERGED_GENCOVE_AXIOM_QC6.fam", header = F,  col.names =c("FID","IID","F","M","sex","size"))
stuck <- read.table("STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.fam", header = F,  col.names =c("FID","IID","F","M","sex","size"))
ccd <- read.table("ALLFAM_DA_MERGED_GENCOVE_AXIOM_QC6.fam", header = F,  col.names =c("FID","IID","F","M","sex","size"))

# Combine all IIDs from phenotype files into one vector
all_IIDs <- unique(c(size$IID, ccd$IID, stuck$IID))
famQC5$in_any <- ifelse(famQC5$IID %in% all_IIDs, 1, -9)

#If you want to keep the same format as your original .fam file, just replace the phenotype column with your new one:
famQC5$phe <- famQC5$in_any
famQC5$in_any <- NULL

write.table(famQC5, "CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC5.fam",quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)

##make the allcontrolFAM:
plink --dog --make-bed --allow-no-sex --prune --geno 0.05 --mind 0.05 --maf 0.01 --hwe 'midp' 0.00000000000000000001 --bed DA_MERGED_GENCOVE_AXIOM_QC5.bed --bim DA_MERGED_GENCOVE_AXIOM_QC5.bim --fam CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC5.fam --out CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6 
#9876571 variants and 2737 dogs pass filters and QC.

#plink files to use: /proj/snic2022-6-229/private/Darwin_Affy/NEW_DA_IMP_MERGED/2024-10-11/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6

#making input plink files for all chr:
for chr in {1..38}; do
    plink --dog \
        --bfile CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6 \
        --chr $chr \
        --make-bed \
        --out CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.$chr
done

###################################################################