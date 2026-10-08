###############################################################################
#####  GENCOVE & AXIOM  imputed files - QC before merging:                ##### 
###############################################################################

#GENCOVE IMP - use output from 01_mapping/3_LowPass_imputation/07_Merge_chr_plink.sh -> DA_IMP_GENCOVE_QC
plink --dog --keep-allele-order --geno 0.05 --mind 0.05 --maf 0.001 --hwe 'midp' 0.00000000000000000001 --bfile DA_IMP_GENCOVE_ALLCHR --make-bed --out DA_IMP_GENCOVE_QC
#--hwe: 61417 variants removed due to Hardy-Weinberg exact test. 1e-20
#854063 variants removed due to minor allele threshold(s)
#13358649 variants and 3044 dogs pass filters and QC.

#AXIOM IMP - use output from 02_data_processing/2_Axiom_imputation/03_axiomimpt_to_plink.sh -> DA_AFFYimp_ALLCHR
plink --allow-no-sex  --dog --geno 0.05 --keep-allele-order --mind 0.05 --maf 0.001 --bfile DA_AFFYimp_ALLCHR --make-bed --out DA_AFFY_IMP --hwe 'midp' 0.00000000000000000001
#--hwe: 116 variants removed due to Hardy-Weinberg exact test. 1e-20
#14042321 variants and 411 dogs pass filters and QC.

### put files in path/to/ and work in that folder...

###############################################################################
#####   Merging    GENCOVE & AXIOM  2024-02                               ##### 
###############################################################################


#merge-mode 6 = (no merge) Report all mismatching calls.
plink --dog --keep-allele-order --bfile DA_IMP_GENCOVE_QC  --bmerge DA_AFFY_IMP --make-bed --merge-mode 6 --merge-equal-pos --out DA_MERGED_GENCOVE_AXIOM
##use --merge-equal-pos. (This will fail if any of the same-position variant pairs do not have matching allele names.) 
#nothing in diff file.

#merge-mode 4 = never overwrite
#merge-mode 1 = (default) Ignore missing calls, otherwise set mismatches to missing not sure if I should run this or merge-mode 4? I run merge-mode=1
plink --dog --keep-allele-order --bfile DA_IMP_GENCOVE_QC  --bmerge DA_AFFY_IMP --make-bed --merge-mode 1 --merge-equal-pos --out DA_MERGED_GENCOVE_AXIOM


#qc merged output: 
plink --dog --bfile DA_MERGED_GENCOVE_AXIOM --make-bed  --geno 0.05 --maf 0.001  --hwe 'midp' 0.00000000000000000001 --out DA_MERGED_GENCOVE_AXIOM_QC
#2275000 variants removed due to missing genotype data (--geno).
#--hwe: 37542 variants removed due to Hardy-Weinberg exact test.
#0 variants removed due to minor allele threshold(s)
#12525443 variants and 3455 dogs pass filters and QC.


#run flipscan on merged qc:ed output:
plink --dog --bfile DA_MERGED_GENCOVE_AXIOM_QC --allow-no-sex --make-pheno DA_IMP_GENCOVE.fam '*' --flip-scan --out flipped
#--flip-scan: 7163 variants with at least one negative LD match.

awk '$10 != "NA"' flipped.flipscan > flipped.flipscan_noNA #7164
awk '{print $2}' flipped.flipscan_noNA > exclude_SNPs.txt

#check maf for flipped markers in both sets :
plink --dog --bfile DA_IMP_GENCOVE_QC --extract exclude_SNPs.txt --freq --out exclude_SNPs_gencove
plink --dog --bfile DA_AFFY_IMP --extract exclude_SNPs.txt --freq --out exclude_SNPs_axiom

#just exclude these SNPs before in both datasets before merging again:
plink --dog --bfile DA_IMP_GENCOVE_QC --exclude exclude_SNPs.txt --make-bed --out DA_IMP_GENCOVE_QC_tmp
plink --dog --bfile DA_AFFY_IMP --exclude exclude_SNPs.txt --make-bed --out DA_AFFY_IMP_tmp

#merge again: use merge-mode 1 and merge-equal positions
plink --dog --bfile DA_IMP_GENCOVE_QC_tmp --bmerge DA_AFFY_IMP_tmp --merge-mode 1 --merge-equal-pos --make-bed --out DA_MERGED_GENCOVE_AXIOM_QC2

#flip scan again: flip2
plink --dog --bfile DA_MERGED_GENCOVE_AXIOM_QC2 --allow-no-sex --make-pheno DA_IMP_GENCOVE_QC.fam '*' --flip-scan --out flipped2
awk '$10 != "NA"' flipped2.flipscan > flipped2.flipscan_noNA #264
awk '{print $2}' flipped2.flipscan_noNA > exclude2_SNPs.txt

#flip2: just exclude these SNPs before in both datasets before merging again:
plink --dog --bfile DA_IMP_GENCOVE_QC_tmp --exclude exclude2_SNPs.txt --make-bed --out DA_IMP_GENCOVE_QC_tmp2
plink --dog --bfile DA_AFFY_IMP_tmp --exclude exclude2_SNPs.txt --make-bed --out DA_AFFY_IMP_tmp2

#merge again: use merge-mode 1 and merge-equal positions
plink --dog --bfile DA_IMP_GENCOVE_QC_tmp2 --bmerge DA_AFFY_IMP_tmp2 --merge-mode 1 --merge-equal-pos --make-bed --out DA_MERGED_GENCOVE_AXIOM_QC3

#flip scan again: flip3
plink --dog --bfile DA_MERGED_GENCOVE_AXIOM_QC3 --allow-no-sex --make-pheno DA_IMP_GENCOVE_QC.fam '*' --flip-scan --out flipped3
#--flip-scan: 0 variants with at least one negative LD match.



##########################################################################################
######## final merged dataset - do qc #####
##########################################################################################


# filter out SNP with low freq and low genotyping rate hwe-1e-20-midp-keep-fewhet(check this later: --ld <variant ID> <variant ID> ['dosage'] ['hwe-midp'])
plink --dog --maf 0.001 --allow-no-sex --geno 0.05 --hwe 0.00000000000000000001 "midp"  --bfile DA_MERGED_GENCOVE_AXIOM_QC3 --make-bed --out DA_MERGED_GENCOVE_AXIOM_QC4
#2275000 variants removed due to missing genotype data (--geno).
#--hwe: 37481 variants removed due to Hardy-Weinberg exact test
#12518078 variants and 3455 dogs pass filters and QC.



******************* start R **************************

# LOAD LIBRARIES
require(tidyverse)
library(dplyr)
setwd("/path/to/")

##ready input files are provided in 02_data_processing/data!

####load QC:ed and modified cf4_fam file with updated dogIDs####
#Here, fam with both sampleID and dogID:
fam <- read.csv("DA_MERGED_GENCOVE_AXIOM_QC3modi.fam", sep= "\t", header=T)

#create list of mean depths from mapping (01_mapping/1_High-LowPass_mapping) and extract_depth_autosomes_and_x.sh
DA_all <- read.csv("DarwinsDogs_N-3285_bam_meandepths.txt", sep= "\t", header=F, col.names = c("sampleID","meanDepthALL"))
DA_all$IID <- as.character(DA_all$sampleID)
DA_merged <- left_join(fam, DA_all, by ="IID")
nrow(subset(DA_merged, is.na(DA_merged$meanDepthALL))) #454 NA -- axiom mainly

#Merge with known duplicates
DAdup <- read.csv("Duplicates_to_keep_DarwinsDogs_N-9_21samples_bam_meandepths.txt", sep= "\t",header=F, col.names = c("Sample_run_id","depth",	"sampleID",	"keep_or_exclude","IID"	,"dogID",	"extra_column"))
DA_merged2 <- left_join(DA_merged, DAdup, by ="IID")

subset(DA_merged2, !is.na(DA_merged2$keep_or_exclude))
nrow(DA_merged2)#3455

DA_merged3 <- DA_merged2 %>%
  mutate(DEPTH_filter=ifelse(meanDepthALL<0.3, 0, 1))
head(DA_merged3)
table(DA_merged3$DEPTH_filter, exclude=NULL)

DA_merged4 <- DA_merged3 %>%
  mutate(DEPTH_filter2=ifelse(is.na(DEPTH_filter), 1, DEPTH_filter))

head(DA_merged4)
table(DA_merged4$DEPTH_filter2, exclude=NULL)#the 0:s are the ones with depth below 0.3
DA_merged5 <- subset(DA_merged4, DA_merged4$DEPTH_filter2=="1")
nrow(DA_merged5)

DA_merged6 <- DA_merged5 %>%
  mutate(keep_or_exclude2=ifelse(is.na(keep_or_exclude), 1, keep_or_exclude))
table(DA_merged6$keep_or_exclude2, exclude=NULL)
DA_merged7 <- subset(DA_merged6, !DA_merged6$keep_or_exclude2=="exclude")
nrow(DA_merged7)
table(DA_merged7$keep_or_exclude2, exclude=NULL)
summary(DA_merged7$meanDepthALL, exclude=NULL)# the keep/duplicates have depth in another column. Checked them manually.

head(DA_merged7)
DA_merged7$qc_pass <- 1

####load fam file from QC4 to merge with the info from QC:ed fam based on depth and duplicates####
famQC4 <- read.csv("DA_MERGED_GENCOVE_AXIOM_QC4.fam", sep= " ", header=F, col.names =c("FID","IID","F","M","sex","phe"))
head(famQC4)

QC4updated <- left_join(famQC4, DA_merged7, by ="FID")
head(QC4updated)
table(QC4updated$qc_pass, exclude=NULL)

#give the ones without a qc_pass =1 (i.e. Okay) a 2 to know which to exclude.
QC4updated2 <- QC4updated %>%
  mutate(qc_pass2=ifelse(is.na(qc_pass), 2, qc_pass))

#make a list of dogs to exclude:
famQC4_exclude <- subset(QC4updated2, QC4updated2$qc_pass2=="2")
nrow(famQC4_exclude)
famQC4_exclude$FID

list_exclude <- data.frame(FID= famQC4_exclude$FID, IID = famQC4_exclude$IID.x) 
#write.table(list_exclude, file= "list_exclude.txt", row.names =F, sep = "\t", col.names= F, quote = F )

******************* end R **************************

#use list to exclude bases on low depth (below 0.3) and known duplicates: 
plink --dog --allow-no-sex --remove list_exclude.txt --bfile DA_MERGED_GENCOVE_AXIOM_QC4 --make-bed --out DA_MERGED_GENCOVE_AXIOM_QC4B
#wc -l DA_MERGED_GENCOVE_AXIOM_QC4B.fam 3341

#The samples are now removed in DA_MERGED_GENCOVE_AXIOM_QC4B

# pruning to run KING:
plink --dog --bfile DA_MERGED_GENCOVE_AXIOM_QC4B --indep-pairwise 50 5 0.2 --out pruned
#Pruning complete.  10113159 of 12518078 variants removed.

plink --dog --pca --bfile DA_MERGED_GENCOVE_AXIOM_QC4B --allow-no-sex --extract pruned.prune.in --out prunedDA_MERGED_GENCOVE_AXIOM_QC4B --make-bed 
#3341 dogs pass filters and QC

#run KING: 
module load KING

king -b prunedDA_MERGED_GENCOVE_AXIOM_QC4B.bed --related --rplot --prefix prunedDA_MERGED_GENCOVE_AXIOM_QC4B --sexchr 39  
#Stages 1&2 (with 32768 SNPs): 144 pairs of relatives are detected (with kinship > 0.1250)
                               Screening ends at Fri Oct 11 15:29:19 2024
  Final Stage (with 2404919 SNPs): 64 pairs of relatives (up to 1st-degree) are confirmed
                               Inference ends at Fri Oct 11 15:29:19 2024

Relationship summary (total relatives: 0 by pedigree, 78 by inference)
        	MZ	PO	FS	2nd	3rd	4th
  =========================================================
  Inference	12	10	42	12	2	0
  


FID1	ID1	DogID	Depth	FID2	ID2	DogID	Part of AXIOM-GENCOVE CONCORDANCE CHECK	Remove Axiom	Depth	N_SNP	HetHet	IBS0	HetConc	HomIBS0	Kinship	IBD1Seg	IBD2Seg	PropIBD	InfType
31001510014690	31001510014690	10567	4.19842	SRR13235184	SRR13235184	10567b	Yes	31001510014690		2404919	0.0627	0.0001	0.9553	0.0031	0.4873	0.0180	0.9585	0.9675	Dup/MZ
31001602030662	31001602030662	10915	1.09361	SRR13233009	SRR13233009	10915b	Yes	31001602030662		2404919	0.0832	0.0000	0.9467	0.0015	0.4860	0.0419	0.9261	0.9471	Dup/MZ
31001602031870	31001602031870	2130	1.37274	SRR13233881	SRR13233881	2130b	Yes	31001602031870		2404919	0.0644	0.0001	0.9578	0.0030	0.4881	0.0240	0.9432	0.9552	Dup/MZ
31001602110587	31001602110587	12129	2.36264	SRR13235060	SRR13235060	12129b	Yes	31001602110587		2404919	0.0770	0.0000	0.9468	0.0022	0.4857	0.0414	0.9321	0.9528	Dup/MZ
31001602111445	31001602111445	7603	0.883885	SRR13233261	SRR13233261	7603b	Yes	31001602111445		2404919	0.0778	0.0001	0.9013	0.0044	0.4729	0.2206	0.7242	0.8345	Dup/MZ
31019073905070	31019073905070	39690	0.793124	31211112308186	31211112308186	30695				2404919	0.0557	0.0000	0.9357	0.0003	0.4831	0.0910	0.8866	0.9321	Dup/MZ
SRR13233291	SRR13233291	157	N/A	SRR13233722	SRR13233722	157b				2404919	0.0893	0.0000	0.9894	0.0000	0.4973	0.0000	0.9963	0.9963	Dup/MZ
SRR13233291	SRR13233291	157	N/A	SRR13235156	SRR13235156	157c				2404919	0.0881	0.0000	0.9737	0.0003	0.4932	0.0161	0.9572	0.9653	Dup/MZ
SRR13233665	SRR13233665	158b	N/A	SRR13233711	SRR13233711	158d				2404919	0.0830	0.0000	0.9642	0.0004	0.4907	0.0100	0.9681	0.9731	Dup/MZ
SRR13233665	SRR13233665	158b	N/A	SRR13235157	SRR13235157	158e				2404919	0.0812	0.0000	0.9427	0.0007	0.4850	0.0453	0.9075	0.9301	Dup/MZ
SRR13233711	SRR13233711	158d	N/A	SRR13235157	SRR13235157	158e				2404919	0.0823	0.0000	0.9567	0.0004	0.4887	0.0312	0.9310	0.9466	Dup/MZ
SRR13233722	SRR13233722	157b	N/A	SRR13235156	SRR13235156	157c				2404919	0.0880	0.0000	0.9718	0.0004	0.4927	0.0184	0.9512	0.9604	Dup/MZ


####
The first five pairs are axiom - genocove duplicates, part of concordance check (five pairs)
Sixth row have different dog ID, maybe monozygotic or just duplicates, remove both. Cannot find 
Dog 157 and 158 - check these, cannot find among depth file. These are found among the 1844 lowpass, keeping the ones with largest bytes. 

#remove second round of duplicates: 11 duplicates
plink --dog --allow-no-sex --remove list_exclude_secondexclusion.txt --bfile DA_MERGED_GENCOVE_AXIOM_QC4B --make-bed --out DA_MERGED_GENCOVE_AXIOM_QC4C

####double check with king:
module load KING
#use same pruned snp set
plink --dog --pca --bfile DA_MERGED_GENCOVE_AXIOM_QC4C --allow-no-sex --extract pruned.prune.in --out prunedDA_MERGED_GENCOVE_AXIOM_QC4C --make-bed 
king -b prunedDA_MERGED_GENCOVE_AXIOM_QC4C.bed --related --rplot --prefix prunedDA_MERGED_GENCOVE_AXIOM_QC4C --sexchr 39  

Relationship summary (total relatives: 0 by pedigree, 63 by inference)
        	MZ	PO	FS	2nd	3rd	4th
  =========================================================
  Inference	0	9	41	11	2	0


############################################
#change fam file from sampleID to dogID!!
#remove duplicate dog 3094, not found among duplicates. Unsure what dog it is. remove:
sampleID	dogID
3094		3094b
SRR13233626	3094

plink --dog --bed DA_MERGED_GENCOVE_AXIOM_QC4C.bed --bim DA_MERGED_GENCOVE_AXIOM_QC4C.bim --fam modi_DA_MERGED_GENCOVE_AXIOM_QC4C_forDogIDchange.fam --remove list_exclude_thirdexclusion.txt  --make-bed --out DA_MERGED_GENCOVE_AXIOM_QC5
#12518078 variants and 3328 dogs pass filters and QC.
#created DA_MERGED_GENCOVE_AXIOM_QC5

******************* start R **************************

###load QC:ed and modified cf4_fam file with updated dogIDs and duplicates and low coverage dogs removed!!!!!!####
#2024-10-14:#in total dogs removed: 106 samples with depth below 0.3, 8 + 11 samples removed due to duplicated samples, 
#and two samples with the same dogID (3094 that did look like a duplicate)
famQC5 <- read.csv("DA_MERGED_GENCOVE_AXIOM_QC5.fam", sep= " ", header=F, col.names =c("FID","IID","F","M","sex","phe"))
famQC5$IID <- as.character(famQC5$IID)

####load covar file incl. q121 - provided by Vista Feb 8th 2024####
updated_covar <- read.csv("DarwinsDogs_Q121_height_age_sex_batch_GENCOVE_AXIOM_QC4.tsv", sep= "\t", header=T)
updated_covar$IID <- as.character(updated_covar$dog)

data1 <- left_join(famQC5, updated_covar, by="IID")
#write.table(data1, file= "data1.txt", row.names =F, sep = "\t", col.names= F, quote = F )

******************* end R **************************

