

#dogs that are duplicates:
ID1	dogIDnew	batch	ID2	dogIDnew	batch
31001510014690	10567	axiom	SRR13235184	10567b	gencove
31001602030662	10915	axiom	SRR13233009	10915b	gencove
31001602031870	2130	axiom	SRR13233881	2130b	gencove
31001602110587	12129	axiom	SRR13235060	12129b	gencove
31001602111445	7603	axiom	SRR13233261	7603b	gencove


##################################
sbatch_wholedog_genotypeconcordance_v3.sh

#!/bin/bash -l
#SBATCH -A naiss2023-5-396
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 65:00:00
#SBATCH -J genotype_conc_mergeoutput


module load bioinfo-tools GATK bcftools 
module load picard/3.1.1 java/OpenJDK_17+35 

RUN_PATH=path/to/GENOTYPE_CONCORDANCE/WHOLE_DOG
GENCOVE_PATH=path/to/5_LigatedImputedCHRs
AFFY_PATH=path/to/Affy_plink/2024-02-08


#make merged vcf files, lists include qc:ed, per chr, data
cd $GENCOVE_PATH
bcftools concat --file-list DA.chr1-38.vcf.list  -O z -o DA.ALL_CHR.vcf.gz
bcftools index  DA.ALL_CHR.vcf.gz

cd $AFFY_PATH
bcftools concat --file-list DA_AFFY.chr1-38.vcf.list   -O z -o DA_AFFY.ALL_CHR.vcf.gz 
bcftools index DA_AFFY.ALL_CHR.vcf.gz 

#extract shared sites:
cd $RUN_PATH
bcftools isec -e'MAF<0.01' -p $RUN_PATH -n=2 -w1 $GENCOVE_PATH'/DA.ALL_CHR.vcf.gz' $AFFY_PATH'/DA_AFFY.ALL_CHR.vcf.gz' 

bcftools view -T sites.txt $GENCOVE_PATH'/DA.ALL_CHR.vcf.gz' | bcftools sort -Oz -o commonDA.ALL_CHR.vcf.gz
bcftools view -T sites.txt $AFFY_PATH'/DA_AFFY.ALL_CHR.vcf.gz'  | bcftools sort -Oz -o commonDA_AFFY.ALL_CHR.vcf.gz

bcftools index commonDA.ALL_CHR.vcf.gz
bcftools index commonDA_AFFY.ALL_CHR.vcf.gz

#run genotype concordance
 #first pair     
java -jar $PICARD_ROOT/picard.jar GenotypeConcordance -CALL_VCF  commonDA.ALL_CHR.vcf.gz -CALL_SAMPLE SRR13235184 -O wholedog_gc_dogID10567_concordance.vcf -TRUTH_VCF commonDA_AFFY.ALL_CHR.vcf.gz -TRUTH_SAMPLE 31001510014690 

#second pair
java -jar $PICARD_ROOT/picard.jar GenotypeConcordance -CALL_VCF  commonDA.ALL_CHR.vcf.gz -CALL_SAMPLE SRR13233009 -O wholedog_gc_dogID10915_concordance.vcf -TRUTH_VCF commonDA_AFFY.ALL_CHR.vcf.gz -TRUTH_SAMPLE 31001602030662 

#third pair
java -jar $PICARD_ROOT/picard.jar GenotypeConcordance -CALL_VCF  commonDA.ALL_CHR.vcf.gz -CALL_SAMPLE SRR13233881 -O wholedog_gc_dogID2130_concordance.vcf -TRUTH_VCF commonDA_AFFY.ALL_CHR.vcf.gz -TRUTH_SAMPLE 31001602031870 

#fourth pair
java -jar $PICARD_ROOT/picard.jar GenotypeConcordance -CALL_VCF  commonDA.ALL_CHR.vcf.gz -CALL_SAMPLE SRR13235060 -O wholedog_gc_dogID12129_concordance.vcf -TRUTH_VCF commonDA_AFFY.ALL_CHR.vcf.gz -TRUTH_SAMPLE 31001602110587 

#fifth pair
java -jar $PICARD_ROOT/picard.jar GenotypeConcordance -CALL_VCF  commonDA.ALL_CHR.vcf.gz -CALL_SAMPLE SRR13233261 -O wholedog_gc_dogID7603concordance.vcf -TRUTH_VCF commonDA_AFFY.ALL_CHR.vcf.gz -TRUTH_SAMPLE 31001602111445 

#not a pair
java -jar $PICARD_ROOT/picard.jar GenotypeConcordance -CALL_VCF  commonDA.ALL_CHR.vcf.gz -CALL_SAMPLE SRR13235060 -O wholedog_gc_dogNO_PAIR_concordance.vcf -TRUTH_VCF commonDA_AFFY.ALL_CHR.vcf.gz -TRUTH_SAMPLE 31001602111445 

#output files:
##reporting NON_REF_GENOTYPE_CONCORDANCE from:
#wholedog_gc_dogID10567_concordance.vcf.genotype_concordance_contingency_metrics
#wholedog_gc_dogID10915_concordance.vcf.genotype_concordance_contingency_metrics
#wholedog_gc_dogID12129_concordance.vcf.genotype_concordance_contingency_metrics
#wholedog_gc_dogID2130_concordance.vcf.genotype_concordance_contingency_metrics
#wholedog_gc_dogID7603concordance.vcf.genotype_concordance_contingency_metrics
#wholedog_gc_dogNO_PAIR_concordance.vcf.genotype_concordance_contingency_metrics