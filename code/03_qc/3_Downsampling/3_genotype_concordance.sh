#merge the chr.vcf files to one:
bcftools concat /path/to/5_LigatedImputedCHRs/Down_0.1x/Down_0.1x.chr{1..38}.merged.vcf.gz -Oz -o Down_0.1x.ALLCHR.merged.qc_modi.vcf.gz
bcftools concat /path/to/5_LigatedImputedCHRs/Down_0.2x/Down_0.2x.chr{1..38}.merged.vcf.gz -Oz -o Down_0.2x.ALLCHR.merged.qc_modi.vcf.gz
bcftools concat /path/to/5_LigatedImputedCHRs/Down_0.3x/Down_0.3x.chr{1..38}.merged.vcf.gz -Oz -o Down_0.3x.ALLCHR.merged.qc_modi.vcf.gz
bcftools concat /path/to/5_LigatedImputedCHRs/Down_0.4x/Down_0.4x.chr{1..38}.merged.vcf.gz -Oz -o Down_0.4x.ALLCHR.merged.qc_modi.vcf.gz
bcftools concat /path/to/5_LigatedImputedCHRs/Down_0.5x/Down_0.5x.chr{1..38}.merged.vcf.gz -Oz -o Down_0.5x.ALLCHR.merged.qc_modi.vcf.gz
bcftools concat /path/to/5_LigatedImputedCHRs/Down_0.6x/Down_0.6x.chr{1..38}.merged.vcf.gz -Oz -o Down_0.6x.ALLCHR.merged.qc_modi.vcf.gz
bcftools concat /path/to/5_LigatedImputedCHRs/Down_0.7x/Down_0.7x.chr{1..38}.merged.vcf.gz -Oz -o Down_0.7x.ALLCHR.merged.qc_modi.vcf.gz
bcftools concat /path/to/5_LigatedImputedCHRs/Down_0.8x/Down_0.8x.chr{1..38}.merged.vcf.gz -Oz -o Down_0.8x.ALLCHR.merged.qc_modi.vcf.gz
bcftools concat /path/to/5_LigatedImputedCHRs/Down_0.9x/Down_0.9x.chr{1..38}.merged.vcf.gz -Oz -o Down_0.9x.ALLCHR.merged.qc_modi.vcf.gz
bcftools concat /path/to/5_LigatedImputedCHRs/Down_1x/Down_1x.chr{1..38}.merged.vcf.gz -Oz -o Down_1x.ALLCHR.merged.qc_modi.vcf.gz

bcftools index Down_0.1x.ALLCHR.merged.qc_modi.vcf.gz
bcftools index Down_0.2x.ALLCHR.merged.qc_modi.vcf.gz
bcftools index Down_0.3x.ALLCHR.merged.qc_modi.vcf.gz
bcftools index Down_0.4x.ALLCHR.merged.qc_modi.vcf.gz
bcftools index Down_0.5x.ALLCHR.merged.qc_modi.vcf.gz
bcftools index Down_0.6x.ALLCHR.merged.qc_modi.vcf.gz
bcftools index Down_0.7x.ALLCHR.merged.qc_modi.vcf.gz
bcftools index Down_0.8x.ALLCHR.merged.qc_modi.vcf.gz
bcftools index Down_0.9x.ALLCHR.merged.qc_modi.vcf.gz
bcftools index Down_1x.ALLCHR.merged.qc_modi.vcf.gz

#check number of variants:

bcftools +counts 10dogs_passfiltered_SNP.HF.ann.id.vcf.gz

bcftools +counts Down_0.1x.ALLCHR.merged.qc_modi.vcf.gz
bcftools +counts Down_0.2x.ALLCHR.merged.qc_modi.vcf.gz
bcftools +counts Down_0.3x.ALLCHR.merged.qc_modi.vcf.gz
bcftools +counts Down_0.4x.ALLCHR.merged.qc_modi.vcf.gz
bcftools +counts Down_0.5x.ALLCHR.merged.qc_modi.vcf.gz
bcftools +counts Down_0.6x.ALLCHR.merged.qc_modi.vcf.gz
bcftools +counts Down_0.7x.ALLCHR.merged.qc_modi.vcf.gz
bcftools +counts Down_0.8x.ALLCHR.merged.qc_modi.vcf.gz
bcftools +counts Down_0.9x.ALLCHR.merged.qc_modi.vcf.gz
bcftools +counts Down_1x.ALLCHR.merged.qc_modi.vcf.gz

plink --dog --vcf Down_0.1x.ALLCHR.merged.qc_modi.vcf.gz --keep-allele-order --double-id --out Down_0.1x.ALLCHR.merged.qc_modi --make-bed 
#--double-id  causes both family and within-family IDs to be set to the sample ID


/path/to/2B_HighCoverageMapping/3_Joint_calling/GenotypeConcordance_HIGHvsLOW/passfiltered_SNP.HF.ann.id.vcf.gz

bcftools view -S samples_to_extract_down.txt /path/to/2B_HighCoverageMapping/3_Joint_calling/GenotypeConcordance_HIGHvsLOW/passfiltered_SNP.HF.ann.id.vcf.gz -Oz -o 10dogs_passfiltered_SNP.HF.ann.id.vcf.gz
bcftools index 10dogs_passfiltered_SNP.HF.ann.id.vcf.gz

samples_to_extract_down.txt #found in /data/03_qc
SRR13340513
SRR13340515
SRR13340516
SRR13340517
SRR13340519
SRR13340526
SRR13340527
SRR13340530
SRR13340531
SRR13340532


wc -l /path/to/Down_0.1x/sites.txt
wc -l /path/to/Down_0.2x/sites.txt
wc -l /path/to/Down_0.3x/sites.txt
wc -l /path/to/Down_0.4x/sites.txt
wc -l /path/to/Down_0.5x/sites.txt
wc -l /path/to/Down_0.6x/sites.txt
wc -l /path/to/Down_0.7x/sites.txt
wc -l /path/to/Down_0.8x/sites.txt
wc -l /path/to/Down_0.9x/sites.txt
wc -l /path/to/Down_1x/sites.txt

##continue in sbatch:
###############################################
sbatch_genoconcordance_DOWNSAMPLED_0.1x.sh
##############################################

#!/bin/bash -l
#SBATCH -A $proj_id
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 85:00:00
#SBATCH -J gen_concDOWNdogs_down0.1x
#SBATCH --mail-user $user_email
#SBATCH --mail-type=ALL

module load bioinfo-tools bcftools 
module load GATK/4.1.4.1
module load picard/3.1.1 java/OpenJDK_17+35 


RUN_PATH=/path/to/Down_0.1x
DOWN_HIGH_PATH=/path/to/
REFERENCE_PATH=/path/to/GSD1.0_RosY

cd $RUN_PATH

#extract shared sites:
bcftools isec -p $RUN_PATH -n=2 -w1 $DOWN_HIGH_PATH'/10dogs_passfiltered_SNP.HF.ann.id.vcf.gz' $DOWN_HIGH_PATH'/Down_0.1x.ALLCHR.merged.qc_modi.vcf.gz'

bcftools view -T sites.txt $DOWN_HIGH_PATH'/Down_0.1x.ALLCHR.merged.qc_modi.vcf.gz' | bcftools sort -Oz -o COMMON_Down_0.1x.ALLCHR.merged.qc_modi.vcf.gz
bcftools view -T sites.txt $DOWN_HIGH_PATH'/10dogs_passfiltered_SNP.HF.ann.id.vcf.gz' | bcftools sort -Oz -o COMMON_10dogs_passfiltered_SNP.HF.ann.id.vcf.gz

bcftools index COMMON_Down_0.1x.ALLCHR.merged.qc_modi.vcf.gz
bcftools index COMMON_10dogs_passfiltered_SNP.HF.ann.id.vcf.gz

#index using gatk:
gatk IndexFeatureFile -I COMMON_Down_0.1x.ALLCHR.merged.qc_modi.vcf.gz
gatk IndexFeatureFile -I COMMON_10dogs_passfiltered_SNP.HF.ann.id.vcf.gz

#update dictionary since the contig information differs between the vcfs
#https://gatk.broadinstitute.org/hc/en-us/articles/360037425731-UpdateVCFSequenceDictionary
gatk UpdateVCFSequenceDictionary -V COMMON_Down_0.1x.ALLCHR.merged.qc_modi.vcf.gz --replace true  --source-dictionary $REFERENCE_PATH'/UU_Cfam_GSD_1.0_ROSY.dict' --output newcontig_COMMON_Down_0.1x.ALLCHR.merged.qc_modi.vcf.gz 
gatk UpdateVCFSequenceDictionary -V COMMON_10dogs_passfiltered_SNP.HF.ann.id.vcf.gz --replace true  --source-dictionary $REFERENCE_PATH'/UU_Cfam_GSD_1.0_ROSY.dict' --output newcontig_COMMON_10dogs_passfiltered_SNP.HF.ann.id.vcf.gz 

bcftools index newcontig_COMMON_Down_0.1x.ALLCHR.merged.qc_modi.vcf.gz
bcftools index newcontig_COMMON_10dogs_passfiltered_SNP.HF.ann.id.vcf.gz


#run genotype concordance
java -jar $PICARD_ROOT/picard.jar GenotypeConcordance -CALL_VCF  newcontig_COMMON_Down_0.1x.ALLCHR.merged.qc_modi.vcf.gz -CALL_SAMPLE SRR13340513 -O SRR13340513.down0.1x_concordance.vcf -TRUTH_VCF newcontig_COMMON_10dogs_passfiltered_SNP.HF.ann.id.vcf.gz -TRUTH_SAMPLE SRR13340513  
java -jar $PICARD_ROOT/picard.jar GenotypeConcordance -CALL_VCF  newcontig_COMMON_Down_0.1x.ALLCHR.merged.qc_modi.vcf.gz -CALL_SAMPLE SRR13340515 -O SRR13340515.down0.1x_concordance.vcf -TRUTH_VCF newcontig_COMMON_10dogs_passfiltered_SNP.HF.ann.id.vcf.gz -TRUTH_SAMPLE SRR13340515  
java -jar $PICARD_ROOT/picard.jar GenotypeConcordance -CALL_VCF  newcontig_COMMON_Down_0.1x.ALLCHR.merged.qc_modi.vcf.gz -CALL_SAMPLE SRR13340516 -O SRR13340516.down0.1x_concordance.vcf -TRUTH_VCF newcontig_COMMON_10dogs_passfiltered_SNP.HF.ann.id.vcf.gz -TRUTH_SAMPLE SRR13340516  
java -jar $PICARD_ROOT/picard.jar GenotypeConcordance -CALL_VCF  newcontig_COMMON_Down_0.1x.ALLCHR.merged.qc_modi.vcf.gz -CALL_SAMPLE SRR13340517 -O SRR13340517.down0.1x_concordance.vcf -TRUTH_VCF newcontig_COMMON_10dogs_passfiltered_SNP.HF.ann.id.vcf.gz -TRUTH_SAMPLE SRR13340517  
java -jar $PICARD_ROOT/picard.jar GenotypeConcordance -CALL_VCF  newcontig_COMMON_Down_0.1x.ALLCHR.merged.qc_modi.vcf.gz -CALL_SAMPLE SRR13340519 -O SRR13340519.down0.1x_concordance.vcf -TRUTH_VCF newcontig_COMMON_10dogs_passfiltered_SNP.HF.ann.id.vcf.gz -TRUTH_SAMPLE SRR13340519  
java -jar $PICARD_ROOT/picard.jar GenotypeConcordance -CALL_VCF  newcontig_COMMON_Down_0.1x.ALLCHR.merged.qc_modi.vcf.gz -CALL_SAMPLE SRR13340526 -O SRR13340526.down0.1x_concordance.vcf -TRUTH_VCF newcontig_COMMON_10dogs_passfiltered_SNP.HF.ann.id.vcf.gz -TRUTH_SAMPLE SRR13340526  
java -jar $PICARD_ROOT/picard.jar GenotypeConcordance -CALL_VCF  newcontig_COMMON_Down_0.1x.ALLCHR.merged.qc_modi.vcf.gz -CALL_SAMPLE SRR13340527 -O SRR13340527.down0.1x_concordance.vcf -TRUTH_VCF newcontig_COMMON_10dogs_passfiltered_SNP.HF.ann.id.vcf.gz -TRUTH_SAMPLE SRR13340527  
java -jar $PICARD_ROOT/picard.jar GenotypeConcordance -CALL_VCF  newcontig_COMMON_Down_0.1x.ALLCHR.merged.qc_modi.vcf.gz -CALL_SAMPLE SRR13340530 -O SRR13340530.down0.1x_concordance.vcf -TRUTH_VCF newcontig_COMMON_10dogs_passfiltered_SNP.HF.ann.id.vcf.gz -TRUTH_SAMPLE SRR13340530  
java -jar $PICARD_ROOT/picard.jar GenotypeConcordance -CALL_VCF  newcontig_COMMON_Down_0.1x.ALLCHR.merged.qc_modi.vcf.gz -CALL_SAMPLE SRR13340531 -O SRR13340531.down0.1x_concordance.vcf -TRUTH_VCF newcontig_COMMON_10dogs_passfiltered_SNP.HF.ann.id.vcf.gz -TRUTH_SAMPLE SRR13340531  
java -jar $PICARD_ROOT/picard.jar GenotypeConcordance -CALL_VCF  newcontig_COMMON_Down_0.1x.ALLCHR.merged.qc_modi.vcf.gz -CALL_SAMPLE SRR13340532 -O SRR13340532.down0.1x_concordance.vcf -TRUTH_VCF newcontig_COMMON_10dogs_passfiltered_SNP.HF.ann.id.vcf.gz -TRUTH_SAMPLE SRR13340532  


####extract the stats from Genotype concordance:
/path/to/Down_0.1x/SRR13340513.down0.1x_concordance.vcf.genotype_concordance_summary_metrics


#make script: statMove.sh
chmod +x statMove.sh

#!/bin/bash
# To run
#               ./statMove.sh samples_to_extract_down.txt      #
#

while read line1
do

######## Variables to change #########
dog_ID=$line1

# Path to folder to store stat files
outputFolder='/path/to/geno_conc_DOWN'
# Path to input folder
baseFolder02x='/path/to/Down_0.2x'
baseFolder03x='/path/to/Down_0.3x/'
baseFolder04x='/path/to/Down_0.4x/'
baseFolder05x='/path/to/Down_0.5x/'
baseFolder06x='/path/to/Down_0.6x/'
baseFolder07x='/path/to/Down_0.7x/'
baseFolder08x='/path/to/Down_0.8x/'
baseFolder09x='/path/to/Down_0.9x/'
baseFolder1x='/path/to/Down_1x/'

echo "#!/bin/bash -l" > statMove.sc;
echo "#SBATCH -A $proj_id" >> statMove.sc;
echo "#SBATCH -p core -n 1" >> statMove.sc;
echo "#SBATCH -J $dog_ID-statMove" >> statMove.sc;
echo "#SBATCH -t 20:00" >> statMove.sc;

echo "dog_ID=$dog_ID" >> statMove.sc;
echo "outputFolder=$outputFolder" >> statMove.sc;
echo "baseFolder=$baseFolder" >> statMove.sc;

echo "cd $baseFolder02x" >> statMove.sc;
echo "cp $dog_ID.down0.2x_concordance.vcf.genotype_concordance_summary_metrics $outputFolder" >> statMove.sc;

echo "cd $baseFolder03x" >> statMove.sc;
echo "cp $dog_ID.down0.3x_concordance.vcf.genotype_concordance_summary_metrics $outputFolder" >> statMove.sc;

echo "cd $baseFolder04x" >> statMove.sc;
echo "cp $dog_ID.down0.4x_concordance.vcf.genotype_concordance_summary_metrics $outputFolder" >> statMove.sc;

echo "cd $baseFolder05x" >> statMove.sc;
echo "cp $dog_ID.down0.5x_concordance.vcf.genotype_concordance_summary_metrics $outputFolder" >> statMove.sc;

echo "cd $baseFolder06x" >> statMove.sc;
echo "cp $dog_ID.down0.6x_concordance.vcf.genotype_concordance_summary_metrics $outputFolder" >> statMove.sc;

echo "cd $baseFolder07x" >> statMove.sc;
echo "cp $dog_ID.down0.7x_concordance.vcf.genotype_concordance_summary_metrics $outputFolder" >> statMove.sc;

echo "cd $baseFolder08x" >> statMove.sc;
echo "cp $dog_ID.down0.8x_concordance.vcf.genotype_concordance_summary_metrics $outputFolder" >> statMove.sc;

echo "cd $baseFolder09x" >> statMove.sc;
echo "cp $dog_ID.down0.9x_concordance.vcf.genotype_concordance_summary_metrics $outputFolder" >> statMove.sc;

echo "cd $baseFolder1x" >> statMove.sc;
echo "cp $dog_ID.down1x_concordance.vcf.genotype_concordance_summary_metrics $outputFolder" >> statMove.sc;

echo "cd $baseFolderALLx" >> statMove.sc;
echo "cp $dog_ID.downALLx_concordance.vcf.genotype_concordance_summary_metrics $outputFolder" >> statMove.sc;



sbatch statMove.sc;

done<samples_to_extract_down.txt 

###########################

chmod +x extract_cov.sh

./extract_cov.sh

#!/bin/bash

output_file="combined_concordance.txt" #found in /data/03_qc

# Clear the output file if it exists
> $output_file

# Use a pattern that matches the specific filename format
for file in SRR*.genotype_concordance_summary_metrics; do
  echo "Processing file: $file"
  
  # Extract sample ID and coverage from filename
  sample_info=$(echo "$file" | grep -oP 'SRR\d+\.down\d+\.\d+x')
  
  # Print the first few lines of the file (for debugging)
  echo "File contents:"
  head -n 20 "$file"
  
  # Extract value from the SNP row (adjust this based on the actual file content)
  value=$(awk '!/^#/ && $1=="SNP" {print $13; exit}' "$file")
  
  echo "Extracted value: $value"
  
  # Combine sample info with extracted value and append to output file
  echo -e "${sample_info}\t${value}" >> $output_file
  
  echo "--------------------------------"
done

# Check if any files were processed
if [ ! -s "$output_file" ]; then
  echo "No matching files found or no data extracted."
  exit 1
fi

echo "Processing complete. Results saved in $output_file"