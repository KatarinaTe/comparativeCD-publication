####
#start with data1 from gencove_axiom_merging_filtering.R

data1 <- read.table("data1.txt",header= T,   sep= "\t")
head(data1)

####dog CD factor scores - produced in the script (EFA and IRT): 1_Survey_data/Darwin_dog-CD_3-factorsolution.R ####
F1_CCD3F <- read.table("F1_CCD3F_240213.txt",header= T,   sep= "\t")
F1_CCD3F$IID <- as.character(F1_CCD3F$dog)

F2_CCD3F <- read.table("F2_CCD3F_240213.txt",header= T,   sep= "\t")
F2_CCD3F$IID <- as.character(F2_CCD3F$dog)

F3_CCD3F <- read.table("F3_CCD3F_240213.txt",header= T,   sep= "\t")
F3_CCD3F$IID <- as.character(F3_CCD3F$dog)


####CCDF1 - merge fam with factor scores####

####merge all factors and items - make: 14CCD_items_phenofile####
CCDitems <- left_join(data1, F1_CCD3F, by="IID")
CCDitems1 <- left_join(CCDitems, F2_CCD3F, by="IID")
data5 <- left_join(CCDitems1, F3_CCD3F, by="IID")

data6 <- data.frame(IID = data5$IID, age = data5$age.x, sex = data5$sex_numeric.x, F1=data5$F1,F1_SE=data5$F1_SE, F2=data5$F2, 
                    F2_SE=data5$F2_SE, F3=data5$F3, F3_SE=data5$F3_SE, 
                    item7=data5$item7, item153=data5$item153, item154=data5$item154,
                    item155=data5$item155, item93=data5$item93, item95=data5$item95, item150=data5$item150,
                    item145=data5$item145, item146=data5$item146, item147=data5$item147, item148=data5$item148, 
                    item149=data5$item149, item151=data5$item151, item152=data5$item152) 


#write.table(data6, file= "data6_2024-10-14.txt", row.names =F, sep = "\t", col.names= T, quote = F )
data6 <- read.csv("data6_2024-10-14.txt", sep= "\t", header=T)
nrow(data6)#3328

####get the allfam dogs with age and sex####
allfam <- read.csv("fam_2584dogs.txt", sep= " ", header=T)
allfam$IID <- as.character(allfam$IID)

all_fam_sex_age <- left_join(allfam, data6, by="IID")
#write.table(all_fam_sex_age, file= "ALLFAM_sex_age_phenos.txt", row.names =F, sep = "\t", col.names= T, quote = F )
#all_fam_sex_age <- read.csv("all_fam_sex_age.txt", sep= "\t", header=T)

####CCDF1 - create dataset####
summary(data6$F1, exclude=NULL) #NA=342
summary(data6$F1_SE, exclude=NULL)

#filter on NA for Factor score

#filter on itemNA count:
head(data6)
data6$na_count_F1 <- apply(data6[,10:13], 1, function(data6) {sum(is.na(data6))})
plot(data6$na_count_F1, xlab="N",pch=19, ylab="NA count F1", title("CCDF1"))#only plots the ones with a factor score

CCDF1b <- data6 %>%
  mutate(F1modi=ifelse(na_count>=3,NA,F1))
table(CCDF1b$na_count, exclude=NULL)
summary(CCDF1b$F1)
summary(CCDF1b$F1modi)

#filter on age-NA:
CCDF1c <- CCDF1b %>%
  mutate(F1modi2=ifelse(is.na(age),NA,F1modi))
table(CCDF1c$F1modi2, exclude=NULL) #NA: 899
nrow(CCDF1c)

#prepare for fam:
CCDF1d <- CCDF1c %>%
  mutate(F1modi3=ifelse(is.na(F1modi2),-9,F1modi2))
head(CCDF1d)
summary(CCDF1d$F1modi3, exclude=NULL) 
table(CCDF1d$sex,CCDF1c$F1modi2,  exclude=NULL) 
CCDF1d$M <- 0

#make fam from CCDF1:
CCDF1_gencove_axiom_fam <- data.frame(FID= CCDF1d$IID, IID = CCDF1d$IID, F=CCDF1d$M, M=CCDF1d$M, 
                                      sex=CCDF1d$sex, phe=CCDF1d$F1modi3) 
head(CCDF1_gencove_axiom_fam)
nrow(CCDF1_gencove_axiom_fam)
#write.table(CCDF1_gencove_axiom_fam, file= "CCDF1_DA_MERGED_GENCOVE_AXIOM_QC5.fam", row.names =F, sep = "\t", col.names= F, quote = F )


#########################################
#TERMINAL run plink on uppmax and produce pruned set QC6:
#make CCDF1 qc:e set QC6:
plink --dog --make-bed --prune --maf 0.01 --bed DA_MERGED_GENCOVE_AXIOM_QC5.bed --bim DA_MERGED_GENCOVE_AXIOM_QC5.bim --fam CCDF1_DA_MERGED_GENCOVE_AXIOM_QC5.fam --out CCDF1_DA_MERGED_GENCOVE_AXIOM_QC6 
#9848182 variants and 2509 dogs pass filters and QC.

awk '{print $1,$2,$6}'  CCDF1_DA_MERGED_GENCOVE_AXIOM_QC6.fam > CCDF1_DA_MERGED_GENCOVE_AXIOM_QC6.phen
#######################################

famQC6 <- read.csv("CCDF1_DA_MERGED_GENCOVE_AXIOM_QC6.fam", sep= " ", header=F, col.names=c("FID","IID","F","M","sex","pheno"))
nrow(famQC6)
famQC6$IID <- as.character(famQC6$IID)
famQC6age <- left_join(famQC6, data4,  by="IID")
nrow(famQC6age)
table(famQC6age$sex_numeric, exclude=NULL)

#make qcovar
famQC6age.qcovar <- data.frame(FID= famQC6age$IID, IID = famQC6age$IID, age=famQC6age$age) 
head(famQC6age.qcovar)
table(famQC6age.qcovar$age, exclude=NULL)
#write.table(famQC6age.qcovar, file= "ageCCDF1_DA_MERGED_GENCOVE_AXIOM_QC6.qcovar", row.names =F, sep = "\t", col.names= F, quote = F )

#make covar 
famQC6sex <- famQC6age %>%
  mutate(sex_binary=ifelse(sex_numeric==1,1,0))
head(famQC6sex)
famQC6sex.covar <- data.frame(FID= famQC6sex$IID, IID = famQC6sex$IID, sex=famQC6sex$sex_binary) 
head(famQC6sex.covar)
tail(famQC6sex.covar)
table(famQC6sex.covar$sex, exclude=NULL)
#write.table(famQC6sex.covar, file= "sexCCDF1_DA_MERGED_GENCOVE_AXIOM_QC6.covar", row.names =F, sep = "\t", col.names= F, quote = F )



####CONTINUE with Factor 2  create dataset 2024-10-15#####
summary(data6$F2, exclude=NULL) #NA=457
summary(data6$F2_SE, exclude=NULL)
head(data6)

#filter on NA for Factor score
#filter on itemNA count:
head(data6)
data6$na_count_F2 <- apply(data6[,14:16], 1, function(data6) {sum(is.na(data6))})
plot( data6$na_count_F2, xlab="NA",pch=19, ylab="NA count F2", title("CCDF2"))#only plots the ones with a factor score

CCDF2b <- data6 %>%
  mutate(F2modi=ifelse(na_count_F2>=2,NA,F2))
table(CCDF2b$na_count_F2, exclude=NULL)
summary(CCDF2b$F2, exclude=NULL)
summary(CCDF2b$F2modi, exclude=NULL)

#filter on age-NA:
CCDF2c <- CCDF2b %>%
  mutate(F2modi2=ifelse(is.na(age),NA,F2modi))
table(CCDF2c$F2modi2, exclude=NULL) #NA: 822
nrow(CCDF2c)

#prepare for fam:
CCDF2d <- CCDF2c %>%
  mutate(F2modi3=ifelse(is.na(F2modi2),-9,F2modi2))
head(CCDF2d)
summary(CCDF2d$F2modi3, exclude=NULL) 

table(CCDF2d$sex,CCDF2c$F2modi2,  exclude=NULL) 

CCDF2d$M <- 0

#make fam from CCDF1:
CCDF2_gencove_axiom_fam <- data.frame(FID= CCDF2d$IID, IID = CCDF2d$IID, F=CCDF2d$M, M=CCDF2d$M, 
                                      sex=CCDF2d$sex, phe=CCDF2d$F2modi3) 
head(CCDF2_gencove_axiom_fam)
nrow(CCDF2_gencove_axiom_fam)
#write.table(CCDF2_gencove_axiom_fam, file= "CCDF2_DA_MERGED_GENCOVE_AXIOM_QC5.fam", row.names =F, sep = "\t", col.names= F, quote = F )

############################
#in TERMINAL run plink on uppmax and produce pruned set QC6:
#make qc:ed set QC6:
plink --dog --make-bed --prune --maf 0.01 --bed DA_MERGED_GENCOVE_AXIOM_QC5.bed --bim DA_MERGED_GENCOVE_AXIOM_QC5.bim --fam CCDF2_DA_MERGED_GENCOVE_AXIOM_QC5.fam --out CCDF2_DA_MERGED_GENCOVE_AXIOM_QC6 
#9854593 variants and 2506 dogs pass filters and QC
awk '{print $1,$2,$6}'  CCDF2_DA_MERGED_GENCOVE_AXIOM_QC6.fam > CCDF2_DA_MERGED_GENCOVE_AXIOM_QC6.phen
############################

famQC6 <- read.csv("CCDF2_DA_MERGED_GENCOVE_AXIOM_QC6.fam", sep= " ", header=F, col.names=c("FID","IID","F","M","sex","pheno"))
nrow(famQC6)
famQC6$IID <- as.character(famQC6$IID)
famQC6age <- left_join(famQC6, data4,  by="IID")
nrow(famQC6age)
table(famQC6age$sex_numeric, exclude=NULL)

#make qcovar
famQC6age.qcovar <- data.frame(FID= famQC6age$IID, IID = famQC6age$IID, age=famQC6age$age) 
head(famQC6age.qcovar)
table(famQC6age.qcovar$age, exclude=NULL)
#write.table(famQC6age.qcovar, file= "ageCCDF2_DA_MERGED_GENCOVE_AXIOM_QC6.qcovar", row.names =F, sep = "\t", col.names= F, quote = F )

#make covar 
famQC6sex <- famQC6age %>%
  mutate(sex_binary=ifelse(sex_numeric==1,1,0))
head(famQC6sex)
famQC6sex.covar <- data.frame(FID= famQC6sex$IID, IID = famQC6sex$IID, sex=famQC6sex$sex_binary) 
head(famQC6sex.covar)
tail(famQC6sex.covar)
table(famQC6sex.covar$sex, exclude=NULL)
#write.table(famQC6sex.covar, file= "sexCCDF2_DA_MERGED_GENCOVE_AXIOM_QC6.covar", row.names =F, sep = "\t", col.names= F, quote = F )



####CONTINUE with Factor 3  create dataset 2024-10-15#####
summary(data6$F3, exclude=NULL) #NA=615
summary(data6$F3_SE, exclude=NULL)
head(data6)

#filter on NA for Factor score
#filter on itemNA count:
head(data6)
data6$na_count_F3 <- apply(data6[,17:23], 1, function(data6) {sum(is.na(data6))})
plot( data6$na_count_F3, xlab="N",pch=19, ylab="NA count F3", title("CCDF2"))#only plots the ones with a factor score

CCDF3b <- data6 %>%
  mutate(F3modi=ifelse(na_count_F3>=7,NA,F3))
table(CCDF3b$na_count_F3, exclude=NULL)
summary(CCDF3b$F3, exclude=NULL)
summary(CCDF3b$F3modi, exclude=NULL)

#filter on age-NA:
CCDF3c <- CCDF3b %>%
  mutate(F3modi2=ifelse(is.na(age),NA,F3modi))
table(CCDF3c$F3modi2, exclude=NULL) #NA: 900
nrow(CCDF3c)

#prepare for fam:
CCDF3d <- CCDF3c %>%
  mutate(F3modi3=ifelse(is.na(F3modi2),-9,F3modi2))
head(CCDF3d)
summary(CCDF3d$F3modi3, exclude=NULL) 

table(CCDF3d$sex,  exclude=NULL) #441 without sex

CCDF3d$M <- 0

#make fam from CCDF1:
CCDF3_gencove_axiom_fam <- data.frame(FID= CCDF3d$IID, IID = CCDF3d$IID, F=CCDF3d$M, M=CCDF3d$M, 
                                      sex=CCDF3d$sex, phe=CCDF3d$F3modi3) 
head(CCDF3_gencove_axiom_fam)
nrow(CCDF3_gencove_axiom_fam)
#write.table(CCDF3_gencove_axiom_fam, file= "CCDF3_DA_MERGED_GENCOVE_AXIOM_QC5.fam", row.names =F, sep = "\t", col.names= F, quote = F )

########################################################
#TERMINAL run plink on uppmax and produce pruned set QC6:
#make CCDF3 qc:e set QC6:
plink --dog --make-bed --prune --maf 0.01 --bed DA_MERGED_GENCOVE_AXIOM_QC5.bed --bim DA_MERGED_GENCOVE_AXIOM_QC5.bim --fam CCDF3_DA_MERGED_GENCOVE_AXIOM_QC5.fam --out CCDF3_DA_MERGED_GENCOVE_AXIOM_QC6 
#9874934 variants and 2428 dogs pass filters and QC.
awk '{print $1,$2,$6}'  CCDF3_DA_MERGED_GENCOVE_AXIOM_QC6.fam > CCDF3_DA_MERGED_GENCOVE_AXIOM_QC6.phen
########################################################

famQC6 <- read.csv("CCDF3_DA_MERGED_GENCOVE_AXIOM_QC6.fam", sep= " ", header=F, col.names=c("FID","IID","F","M","sex","pheno"))
nrow(famQC6) #2428
famQC6$IID <- as.character(famQC6$IID)
famQC6age <- left_join(famQC6, data4,  by="IID")
nrow(famQC6age)
table(famQC6age$sex_numeric, exclude=NULL)

#make qcovar
famQC6age.qcovar <- data.frame(FID= famQC6age$IID, IID = famQC6age$IID, age=famQC6age$age) 
head(famQC6age.qcovar)
table(famQC6age.qcovar$age, exclude=NULL)
#write.table(famQC6age.qcovar, file= "ageCCDF3_DA_MERGED_GENCOVE_AXIOM_QC6.qcovar", row.names =F, sep = "\t", col.names= F, quote = F )

#make covar 
famQC6sex <- famQC6age %>%
  mutate(sex_binary=ifelse(sex_numeric==1,1,0))
head(famQC6sex)
famQC6sex.covar <- data.frame(FID= famQC6sex$IID, IID = famQC6sex$IID, sex=famQC6sex$sex_binary) 
head(famQC6sex.covar)
tail(famQC6sex.covar)
table(famQC6sex.covar$sex, exclude=NULL)
#write.table(famQC6sex.covar, file= "sexCCDF3_DA_MERGED_GENCOVE_AXIOM_QC6.covar", row.names =F, sep = "\t", col.names= F, quote = F )



