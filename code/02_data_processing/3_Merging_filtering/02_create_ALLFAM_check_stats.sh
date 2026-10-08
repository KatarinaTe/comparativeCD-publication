##################################################################
#make an all fam file with all dogs occuring in any of the phenotypes
##################################################################

###ready input files are provided in 02_data_processing/data: fam_2584dogs.txt
#qc:e set QC6:
plink --dog --make-bed --maf 0.01 --bfile DA_MERGED_GENCOVE_AXIOM_QC5 --keep fam_2584dogs.txt --out ALLFAM_DA_MERGED_GENCOVE_AXIOM_QC6 

# pruning
plink --dog --bfile ALLFAM_DA_MERGED_GENCOVE_AXIOM_QC6 --indep-pairwise 50 5 0.2 --out ALLFAM_pruned
#9877574 variants loaded from .bim file.
#2584 dogs (0 males, 0 females, 2584 ambiguous) loaded from .fam.
#Pruning complete.  10113159 of 12518078 variants removed.
# 8771510 of 9877574 variants removed.

##create eigenvectors for PCA plotting (4_Population_structure)
plink --dog --pca --bfile ALLFAM_DA_MERGED_GENCOVE_AXIOM_QC6 --allow-no-sex --extract ALLFAM_pruned.prune.in --out prunedALLFAM_DA_MERGED_GENCOVE_AXIOM_QC6 --make-bed 

#make grm
gcta64 --bfile ALLFAM_DA_MERGED_GENCOVE_AXIOM_QC6 --autosome  --make-grm  --out ALLFAM_DA_MERGED_GENCOVE_AXIOM_QC6 --autosome-num 38

##################################################################
#check SNPs and dogs in each GWAS-dataset:
##################################################################

###terminal

#N dogs:
wc -l CCDF1_DA_MERGED_GENCOVE_AXIOM_QC6.fam
wc -l CCDF2_DA_MERGED_GENCOVE_AXIOM_QC6.fam
wc -l CCDF3_DA_MERGED_GENCOVE_AXIOM_QC6.fam
wc -l CCDitem7_DA_MERGED_GENCOVE_AXIOM_QC6.fam
wc -l CCDitem93_DA_MERGED_GENCOVE_AXIOM_QC6.fam
wc -l CCDitem95_DA_MERGED_GENCOVE_AXIOM_QC6.fam
wc -l CCDitem145_DA_MERGED_GENCOVE_AXIOM_QC6.fam
wc -l CCDitem146_DA_MERGED_GENCOVE_AXIOM_QC6.fam
wc -l CCDitem147_DA_MERGED_GENCOVE_AXIOM_QC6.fam
wc -l CCDitem148_DA_MERGED_GENCOVE_AXIOM_QC6.fam
wc -l CCDitem149_DA_MERGED_GENCOVE_AXIOM_QC6.fam
wc -l CCDitem150_DA_MERGED_GENCOVE_AXIOM_QC6.fam
wc -l CCDitem151_DA_MERGED_GENCOVE_AXIOM_QC6.fam
wc -l CCDitem152_DA_MERGED_GENCOVE_AXIOM_QC6.fam
wc -l CCDitem153_DA_MERGED_GENCOVE_AXIOM_QC6.fam
wc -l CCDitem154_DA_MERGED_GENCOVE_AXIOM_QC6.fam
wc -l CCDitem155_DA_MERGED_GENCOVE_AXIOM_QC6.fam


#N SNPs:
wc -l CCDF1_DA_MERGED_GENCOVE_AXIOM_QC6.bim
wc -l CCDF2_DA_MERGED_GENCOVE_AXIOM_QC6.bim
wc -l CCDF3_DA_MERGED_GENCOVE_AXIOM_QC6.bim
wc -l CCDitem7_DA_MERGED_GENCOVE_AXIOM_QC6.bim
wc -l CCDitem93_DA_MERGED_GENCOVE_AXIOM_QC6.bim
wc -l CCDitem95_DA_MERGED_GENCOVE_AXIOM_QC6.bim
wc -l CCDitem145_DA_MERGED_GENCOVE_AXIOM_QC6.bim
wc -l CCDitem146_DA_MERGED_GENCOVE_AXIOM_QC6.bim
wc -l CCDitem147_DA_MERGED_GENCOVE_AXIOM_QC6.bim
wc -l CCDitem148_DA_MERGED_GENCOVE_AXIOM_QC6.bim
wc -l CCDitem149_DA_MERGED_GENCOVE_AXIOM_QC6.bim
wc -l CCDitem150_DA_MERGED_GENCOVE_AXIOM_QC6.bim
wc -l CCDitem151_DA_MERGED_GENCOVE_AXIOM_QC6.bim
wc -l CCDitem152_DA_MERGED_GENCOVE_AXIOM_QC6.bim
wc -l CCDitem153_DA_MERGED_GENCOVE_AXIOM_QC6.bim
wc -l CCDitem154_DA_MERGED_GENCOVE_AXIOM_QC6.bim
wc -l CCDitem155_DA_MERGED_GENCOVE_AXIOM_QC6.bim

QC6	Dogs	SNPs
CCDF1	2429	9875386
CCDF2	2506	9854593
CCDF3	2428	9874934
CCDitem7	2555	9857137
CCDitem93	2494	9884559
CCDitem95	2498	9885182
CCDitem145	2418	9869715
CCDitem146	2419	9868371
CCDitem147	2414	9866711
CCDitem148	2424	9869166
CCDitem149	2422	9871025
CCDitem150	2417	9867803
CCDitem151	2421	9868912
CCDitem152	2415	9865816
CCDitem153	2417	9868642
CCDitem154	2418	9867972
CCDitem155	2322	9871386



#######################################################
##create stats for the samples used in our analyses: 2026-05-06
setwd("~/Documents/GitHub/comparativeCD/02_data_processing/3_Merging_filtering")
library(ggplot2)
require(tidyverse)
library(dplyr)
library(psych)
library(patchwork)
library(ggtext)

DA_all <- read.csv("~/Documents/GitHub/comparativeCD/02_data_processing/data/DarwinsDogs_N-3285_bam_meandepths.txt", sep= "\t", header=F, col.names = c("sampleID","meanDepthALL"))
head(DA_all)
DA_all$IID <- DA_all$sampleID

summary(DA_all$meanDepthALL, exclude=NULL)
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.000236 0.576519 0.763423 0.807144 0.978977 4.751190 

ALLFAM <- read.csv("~/Documents/GitHub/comparativeCD/04_gwas/data/ALLFAM_DA_MERGED_GENCOVE_AXIOM_QC6.fam", sep= " ", header=F, col.names =c("FID","IID","F","M","sex","phe"))
head(ALLFAM)

depth_allfam <- left_join(ALLFAM, DA_all, by="IID")
summary(depth_allfam$meanDepthALL, exclude=NULL) #373 axiom
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 0.3012  0.5963  0.7764  0.8313  0.9912  4.7512     373 
sd(!is.na(depth_allfam$meanDepthALL))
#[1] 0.3515121

nrow(depth_allfam)
2584-373 #2211

hist_1 <- ggplot(DA_merged, aes(x=meanDepthALL)) +
  geom_histogram(alpha = 1,binwidth=0.01, color="coral3", fill="coral3") +
  scale_fill_manual(values=c("coral3")) +
  scale_color_manual(values=c("coral3")) +
  xlim(-0.1,5)+
  #ylim(0,1800) +
  #scale_x_continuous(breaks = -0.05:0.05, limits = c(-0.05, 0.05)) +
  annotate("text", x = Inf, y = Inf, label = "N=3285", hjust = 1.3, vjust = 2.5, size = 4) +
  labs(title = "ALL DARWIN depth all sites", x = "Depth", y = "Number of dogs")+
  theme(plot.title.position = "plot", 
        plot.title = element_textbox_simple(size = 10, padding = margin(5.5, 5.5, 5.5, 5.5), margin = margin(0, 0, 5.5, 38), 
                                            fill = "cornsilk"),
        axis.title.x = element_textbox_simple(width = NULL, padding = margin(4, 4, 4, 4), margin = margin(4, 0, 0, 0), 
                                              fill = "white"),
        axis.title.y = element_textbox_simple(hjust = 0.65, orientation = "left-rotated", minwidth = unit(1, "in"), maxwidth = unit(2, "in"),
                                              padding = margin(4, 4, 2, 4), margin = margin(0, 0, 2,15), fill = "white"))

hist_4 <- ggplot(depth_allfam, aes(x=meanDepthALL)) +
  geom_histogram(alpha = 1,binwidth=0.01, color="cornflowerblue", fill="cornflowerblue") +
  scale_fill_manual(values=c("cornflowerblue")) +
  scale_color_manual(values=c("cornflowerblue")) +
  xlim(-0.1,5)+
  #ylim(0,1800) +
  #scale_x_continuous(breaks = -0.05:0.05, limits = c(-0.05, 0.05)) +
  annotate("text", x = Inf, y = Inf, label = "N=2211", hjust = 1.3, vjust = 2.5, size = 4) +
  labs(title = "ALLFAM DARWIN depth all sites, mean=0.83, sd=0.35", x = "Depth", y = "Number of dogs")+
  theme(plot.title.position = "plot", 
        plot.title = element_textbox_simple(size = 10, padding = margin(5.5, 5.5, 5.5, 5.5), margin = margin(0, 0, 5.5, 38), 
                                            fill = "cornsilk"),
        axis.title.x = element_textbox_simple(width = NULL, padding = margin(4, 4, 4, 4), margin = margin(4, 0, 0, 0), 
                                              fill = "white"),
        axis.title.y = element_textbox_simple(hjust = 0.65, orientation = "left-rotated", minwidth = unit(1, "in"), maxwidth = unit(2, "in"),
                                              padding = margin(4, 4, 2, 4), margin = margin(0, 0, 2,15), fill = "white"))


hist_1/hist_4


#####