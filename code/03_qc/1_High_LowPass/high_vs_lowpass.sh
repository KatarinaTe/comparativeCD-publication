
###########################################################################################################################
####2024-10-02: redo the genotype concordance - first exclude the variants that did not PASS
###########################################################################################################################

#first exclude the variants that did not PASS
bcftools filter -i 'PASS' filtered_SNP.HF.ann.id.vcf.gz -Oz -o passfiltered_SNP.HF.ann.id.vcf.gz
bcftools view -f 'PASS' filtered_SNP.HF.ann.id.vcf.gz -Oz -o passfiltered_SNP.HF.ann.id.vcf.gz

zcat filtered_SNP.HF.ann.id.vcf.gz | grep -v "#" | wc -l #11010306
zcat passfiltered_SNP.HF.ann.id.vcf.gz | grep -v "#" | wc -l #10499622

bcftools index passfiltered_SNP.HF.ann.id.vcf.gz

###filter out our 11 individuals in 
bcftools view -S samples_to_extract.txt /path/to/5_LigatedImputedCHRs/DA.ALL_CHR.vcf.gz -Oz -o 11dogs_DA.ALL_CHR.vcf.gz
bcftools index 11dogs_DA.ALL_CHR.vcf.gz
#samples_to_extract.txt #available in data/03_qc

#filter on maf in the 11 dogs
bcftools view -e 'MAF<0.01' 11dogs_DA.ALL_CHR.vcf.gz -Oz -o 11dogs_DA.ALL_CHR.maf0.01.vcf.gz
bcftools index 11dogs_DA.ALL_CHR.maf0.01.vcf.gz

bcftools view -e 'MAF<0.01' passfiltered_SNP.HF.ann.id.vcf.gz  -Oz -o passfiltered_SNP.HF.maf0.01.ann.id.vcf.gz 
bcftools index passfiltered_SNP.HF.maf0.01.ann.id.vcf.gz

##continue in sbatch:
sbatch_genoconcordance_11dogs_v2.sh

#!/bin/bash -l
#SBATCH -A naiss2024-5-494
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 85:00:00
#SBATCH -J gen_conc11dogs

module load bioinfo-tools bcftools 
module load GATK/4.1.4.1
module load picard/3.1.1 java/OpenJDK_17+35 

RUN_PATH=path/to/3_Joint_calling/GenotypeConcordance_HIGHvsLOW
REFERENCE_PATH=/path/to/GSD1.0_RosY

cd $RUN_PATH

#extract shared sites:
bcftools isec -p path/to/GenotypeConcordance_HIGHvsLOW -n=2 -w1 passfiltered_SNP.HF.maf0.01.ann.id.vcf.gz 11dogs_DA.ALL_CHR.maf0.01.vcf.gz 
#wc -l sites.txt: 7824389 before (note: when not filted for pass)! Renamed to sites_240910.txt
#Now: 8090162 sites_240910.txt

bcftools view -T sites.txt 11dogs_DA.ALL_CHR.maf0.01.vcf.gz | bcftools sort -Oz -o 11dogs.commonDA_LOW_COV.vcf.gz
bcftools view -T sites.txt passfiltered_SNP.HF.maf0.01.ann.id.vcf.gz | bcftools sort -Oz -o 11dogs.commonDA_HIGH_COV.vcf.gz

bcftools index 11dogs.commonDA_LOW_COV.vcf.gz
bcftools index 11dogs.commonDA_HIGH_COV.vcf.gz

#index using gatk:
gatk IndexFeatureFile -I 11dogs.commonDA_LOW_COV.vcf.gz
gatk IndexFeatureFile -I 11dogs.commonDA_HIGH_COV.vcf.gz

#update dictionary since the contig information differs between the vcfs
#https://gatk.broadinstitute.org/hc/en-us/articles/360037425731-UpdateVCFSequenceDictionary
gatk UpdateVCFSequenceDictionary -V 11dogs.commonDA_LOW_COV.vcf.gz --replace true  --source-dictionary $REFERENCE_PATH'/UU_Cfam_GSD_1.0_ROSY.dict' --output newcontig_11dogs.commonDA_LOW_COV.vcf.gz 
gatk UpdateVCFSequenceDictionary -V 11dogs.commonDA_HIGH_COV.vcf.gz --replace true  --source-dictionary $REFERENCE_PATH'/UU_Cfam_GSD_1.0_ROSY.dict' --output newcontig_11dogs.commonDA_HIGH_COV.vcf.gz 

bcftools index newcontig_11dogs.commonDA_LOW_COV.vcf.gz
bcftools index newcontig_11dogs.commonDA_HIGH_COV.vcf.gz


#run genotype concordance
java -jar $PICARD_ROOT/picard.jar GenotypeConcordance -CALL_VCF  newcontig_11dogs.commonDA_LOW_COV.vcf.gz -CALL_SAMPLE SRR13233711 -O Kaylee.v2_concordance.vcf -TRUTH_VCF newcontig_11dogs.commonDA_HIGH_COV.vcf.gz  -TRUTH_SAMPLE SRR13340530  
java -jar $PICARD_ROOT/picard.jar GenotypeConcordance -CALL_VCF  newcontig_11dogs.commonDA_LOW_COV.vcf.gz -CALL_SAMPLE SRR13233676 -O Kaylee_dup1.v2_concordance.vcf -TRUTH_VCF newcontig_11dogs.commonDA_HIGH_COV.vcf.gz  -TRUTH_SAMPLE SRR13340530  
java -jar $PICARD_ROOT/picard.jar GenotypeConcordance -CALL_VCF  newcontig_11dogs.commonDA_LOW_COV.vcf.gz -CALL_SAMPLE SRR13233654 -O Kaylee_dup2.v2_concordance.vcf -TRUTH_VCF newcontig_11dogs.commonDA_HIGH_COV.vcf.gz  -TRUTH_SAMPLE SRR13340530  
java -jar $PICARD_ROOT/picard.jar GenotypeConcordance -CALL_VCF  newcontig_11dogs.commonDA_LOW_COV.vcf.gz -CALL_SAMPLE SRR13233665 -O Kaylee_dup3.v2_concordance.vcf -TRUTH_VCF newcontig_11dogs.commonDA_HIGH_COV.vcf.gz  -TRUTH_SAMPLE SRR13340530  
java -jar $PICARD_ROOT/picard.jar GenotypeConcordance -CALL_VCF  newcontig_11dogs.commonDA_LOW_COV.vcf.gz -CALL_SAMPLE SRR13233313 -O Clarence.v2_concordance.vcf -TRUTH_VCF newcontig_11dogs.commonDA_HIGH_COV.vcf.gz  -TRUTH_SAMPLE SRR13340532  
java -jar $PICARD_ROOT/picard.jar GenotypeConcordance -CALL_VCF  newcontig_11dogs.commonDA_LOW_COV.vcf.gz -CALL_SAMPLE SRR13233291 -O Beskow.v2_concordance.vcf -TRUTH_VCF newcontig_11dogs.commonDA_HIGH_COV.vcf.gz  -TRUTH_SAMPLE SRR13340531  
java -jar $PICARD_ROOT/picard.jar GenotypeConcordance -CALL_VCF  newcontig_11dogs.commonDA_LOW_COV.vcf.gz -CALL_SAMPLE SRR13233722 -O Beskow_dup1.v2_concordance.vcf -TRUTH_VCF newcontig_11dogs.commonDA_HIGH_COV.vcf.gz  -TRUTH_SAMPLE SRR13340531  
java -jar $PICARD_ROOT/picard.jar GenotypeConcordance -CALL_VCF  newcontig_11dogs.commonDA_LOW_COV.vcf.gz -CALL_SAMPLE SRR13233642 -O Gus.v2_concordance.vcf -TRUTH_VCF newcontig_11dogs.commonDA_HIGH_COV.vcf.gz  -TRUTH_SAMPLE SRR13340529  
java -jar $PICARD_ROOT/picard.jar GenotypeConcordance -CALL_VCF  newcontig_11dogs.commonDA_LOW_COV.vcf.gz -CALL_SAMPLE SRR13233900 -O Esme.v2_concordance.vcf -TRUTH_VCF newcontig_11dogs.commonDA_HIGH_COV.vcf.gz  -TRUTH_SAMPLE SRR13340527  
java -jar $PICARD_ROOT/picard.jar GenotypeConcordance -CALL_VCF  newcontig_11dogs.commonDA_LOW_COV.vcf.gz -CALL_SAMPLE SRR13233889 -O Lucky.v2_concordance.vcf -TRUTH_VCF newcontig_11dogs.commonDA_HIGH_COV.vcf.gz  -TRUTH_SAMPLE SRR13340526  
java -jar $PICARD_ROOT/picard.jar GenotypeConcordance -CALL_VCF  newcontig_11dogs.commonDA_LOW_COV.vcf.gz -CALL_SAMPLE SRR13233643 -O Tui.v2_concordance.vcf -TRUTH_VCF newcontig_11dogs.commonDA_HIGH_COV.vcf.gz  -TRUTH_SAMPLE SRR13340519  
java -jar $PICARD_ROOT/picard.jar GenotypeConcordance -CALL_VCF  newcontig_11dogs.commonDA_LOW_COV.vcf.gz -CALL_SAMPLE SRR13233628 -O Lily.v2_concordance.vcf -TRUTH_VCF newcontig_11dogs.commonDA_HIGH_COV.vcf.gz  -TRUTH_SAMPLE SRR13340517  
java -jar $PICARD_ROOT/picard.jar GenotypeConcordance -CALL_VCF  newcontig_11dogs.commonDA_LOW_COV.vcf.gz -CALL_SAMPLE SRR13233162 -O Finch.v2_concordance.vcf -TRUTH_VCF newcontig_11dogs.commonDA_HIGH_COV.vcf.gz  -TRUTH_SAMPLE SRR13340516  
java -jar $PICARD_ROOT/picard.jar GenotypeConcordance -CALL_VCF  newcontig_11dogs.commonDA_LOW_COV.vcf.gz -CALL_SAMPLE SRR13233942 -O Jamie.v2_concordance.vcf -TRUTH_VCF newcontig_11dogs.commonDA_HIGH_COV.vcf.gz  -TRUTH_SAMPLE SRR13340515  
java -jar $PICARD_ROOT/picard.jar GenotypeConcordance -CALL_VCF  newcontig_11dogs.commonDA_LOW_COV.vcf.gz -CALL_SAMPLE SRR13234818 -O Rudy.v2_concordance.vcf -TRUTH_VCF newcontig_11dogs.commonDA_HIGH_COV.vcf.gz  -TRUTH_SAMPLE SRR13340513  

#output files:
##reporting NON_REF_GENOTYPE_CONCORDANCE from:
#Beskow_dup1.v2_concordance.vcf.genotype_concordance_summary_metrics
#Beskow.v2_concordance.vcf.genotype_concordance_summary_metrics
#Clarence.v2_concordance.vcf.genotype_concordance_summary_metrics
#Esme.v2_concordance.vcf.genotype_concordance_summary_metrics
#Finch.v2_concordance.vcf.genotype_concordance_summary_metrics
#Gus.v2_concordance.vcf.genotype_concordance_summary_metrics
#Jamie.v2_concordance.vcf.genotype_concordance_summary_metrics
#Kaylee_dup1.v2_concordance.vcf.genotype_concordance_summary_metrics
#Kaylee_dup2.v2_concordance.vcf.genotype_concordance_summary_metrics
#Kaylee_dup3.v2_concordance.vcf.genotype_concordance_summary_metrics
#Kaylee.v2_concordance.vcf.genotype_concordance_summary_metrics
#Lily.v2_concordance.vcf.genotype_concordance_summary_metrics
#Lucky.v2_concordance.vcf.genotype_concordance_summary_metrics
#Rudy.v2_concordance.vcf.genotype_concordance_summary_metrics
#Tui.v2_concordance.vcf.genotype_concordance_summary_metrics