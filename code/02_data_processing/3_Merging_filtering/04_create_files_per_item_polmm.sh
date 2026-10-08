###############
*** OBS repeat for each item, ie. 14 items!!!!!!!!!!!!!!!!!!!!!!
###############

module load zstd
module load bioinfo-tools R/4.3.1
module load R_packages/4.3.1

#################
R
require(tidyverse)
library(dplyr)

data6 <- read.csv("data6_2024-10-14.txt", sep= "\t", header=T)
nrow(data6)#3328

####item155 ####
table(data6$item155, exclude=NULL) #NA: 629
#filter on age-NA:
item155 <- data6 %>%
  mutate(item155_modi=ifelse(is.na(age),NA,item155))
table(item155$item155_modi, exclude=NULL) #NA: 913
nrow(item155)
table(item155$item155_modi, item155$sex, exclude=NULL) #441 NA - 
#no need to exclude dogs without sex info, they dont have the phenotype anyway

#prepare for fam:
item155b <- item155 %>%
  mutate(item155_modi2=ifelse(is.na(item155_modi),-9,item155_modi))
head(item155b)
table(item155b$item155_modi2) 
item155b$M <- 0

#make fam from this item:
item155b_fam <- data.frame(FID= item155b$IID, IID = item155b$IID, F=item155b$M, M=item155b$M, 
                           sex=item155b$sex, phe=item155b$item155_modi2) 

write.table(item155b_fam, file= "CCDitem155_DA_MERGED_GENCOVE_AXIOM_QC5.fam", row.names =F, sep = "\t", col.names= F, quote = F )
################################################################################################################################
#in terminal plink:
#make item155 qc:e set QC6:
plink --dog --allow-no-sex --make-bed --prune --maf 0.01 --bed DA_MERGED_GENCOVE_AXIOM_QC5.bed --bim DA_MERGED_GENCOVE_AXIOM_QC5.bim --fam CCDitem155_DA_MERGED_GENCOVE_AXIOM_QC5.fam --out CCDitem155_DA_MERGED_GENCOVE_AXIOM_QC6 

#make a dataset with fewer markers - i.e. LD-pruned.
plink --bfile CCDitem155_DA_MERGED_GENCOVE_AXIOM_QC6 --dog --indep-pairwise 50 5 0.2 --out CCDitem155_indepSNP 
plink --bfile CCDitem155_DA_MERGED_GENCOVE_AXIOM_QC6 --dog --extract CCDitem155_indepSNP.prune.in --make-bed --out CCDitem155_EigenInput

################################################################################################################################################################
################################################################################################################################################################
#in R again: (have two terminal windows open)


#run plink on uppmax and produce pruned set QC6:
famQC6 <- read.csv("CCDitem155_DA_MERGED_GENCOVE_AXIOM_QC6.fam", sep= " ", header=F, col.names=c("FID","IID","F","M","sex","pheno"))
head(famQC6)
table(famQC6$pheno, exclude=NULL)
famQC6$IID <- as.character(famQC6$IID)

#if need to reload:
#data6 <- read.csv("data6_2024-10-14.txt", sep= "\t", header=T) 
#head(data6)

famQC6b <- left_join(famQC6, data6,  by="IID")
head(famQC6b)

famQC6c <- famQC6b %>%
  mutate(sex_binary=ifelse(sex.x==1,1,0)) #works because no sex NA!!!

QC6grab <- data.frame(IID = famQC6c$IID, item155=famQC6c$item155, sex=famQC6c$sex_binary, age=famQC6c$age) 
head(QC6grab)
nrow(QC6grab) #2322
write.table(QC6grab, file= "grabCCDitem155_DA_MERGED_GENCOVE_AXIOM_QC6.txt", row.names =F, sep = "\t", col.names= T, quote = F )

################################################################################################################################################################
