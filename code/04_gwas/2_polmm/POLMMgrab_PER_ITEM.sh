
#do this in R on uppmax. No need to transfer.
module load zstd
module load bioinfo-tools R/4.3.1
module load R_packages/4.3.1


#OBS!!! #Here example for item155: Change and run for each item:

########################################################################################################################
#Rscript polmm_CCDitem155.R ##############################
R
library(GRAB)
setwd("/path/to/")


PhenoData <- read.table("grabCCDitem155_DA_MERGED_GENCOVE_AXIOM_QC6.txt", header = T, sep = "") 
PhenoData$item155_factor <- as.factor(PhenoData$item155)
GenoFile = "CCDitem155_DA_MERGED_GENCOVE_AXIOM_QC6.bed"
prunedGeno = "CCDitem155_EigenInput.bed" #ca 1 million SNPs

# Step 1: fit a null model - using pruned genodata to make grm
obj.POLMM_item155_FULLGRM = GRAB.NullModel(formula = item155_factor ~ sex + age,
                           data = PhenoData, 
                           subjData = PhenoData$IID, 
                           method = "POLMM", 
                           traitType = "ordinal",
                           GenoFile = prunedGeno,
                           control = list(showInfo = FALSE, 
                                          LOCO = FALSE, 
                                          tolTau = 0.2, 
                                          tolBeta = 0.1))                       


#step 2 full genofile:
OutputDir = "./"
OutputFile = paste0(OutputDir, "/simuMarkerOutput_POLMM_item155_FULLGRM_fullGeno.txt")

GRAB.Marker(obj.POLMM_item155_FULLGRM, GenoFile = GenoFile, OutputFile = OutputFile)


########################################################################################################################
#polmm sumstats was modified like this:
awk -F'\t' ' NR==1 {print "CHR\tPosition\tChr\tbp\tA1\tA2\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8; next} { split($1,a,":"); split($2,b,":");print a[1] "\t" a[2] "\t" b[1] "\t" b[2] "\t" b[3] "\t" b[4] "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8}' simuMarkerOutput_POLMM_item155_FULLGRM_fullGeno.txt > modi_simuMarkerOutput_POLMM_item155_FULLGRM_fullGeno.txt
########################################################################################################################
