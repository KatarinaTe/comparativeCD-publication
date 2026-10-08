

##use input paths
split -l 1 DA_downsampling_bam_imput.txt down 

## Step 1. mapping to make bams. OBS new EDITS  (start with BAM not FASTQ)

for i in down*;do echo $i; sbatch dog10k_mapping_start_at_bams.sh $i; done	#obs updated dog10k_mapping script, see 03_qc/data


###Step 2. Genotype likelihood
# Run genotypelikelihood: in batches for each downsample-batch do as regular mapping, update path to bams and run script: ./02_GenotypeLikelihood-Down0.1x.sh CHR1-38.txt etc.
--bam-list ${ScriptsFolder}/DownSampledBams.0.1x.txt
--bam-list ${ScriptsFolder}/DownSampledBams.0.2x.txt
--bam-list ${ScriptsFolder}/DownSampledBams.0.3x.txt
--bam-list ${ScriptsFolder}/DownSampledBams.0.4x.txt
--bam-list ${ScriptsFolder}/DownSampledBams.0.5x.txt
--bam-list ${ScriptsFolder}/DownSampledBams.0.6x.txt
--bam-list ${ScriptsFolder}/DownSampledBams.0.7x.txt
--bam-list ${ScriptsFolder}/DownSampledBams.0.8x.txt
--bam-list ${ScriptsFolder}/DownSampledBams.0.9x.txt
--bam-list ${ScriptsFolder}/DownSampledBams.1x.txt

#e.g. ${ScriptsFolder}/DownSampledBams.0.1x.txt
#/path/to/SRR13340513.0.1x/SRR13340513.0.1x.sorted.merged.MarkDups.BQSR.bam
#/path/to/SRR13340515.0.1x/SRR13340515.0.1x.sorted.merged.MarkDups.BQSR.bam
#/path/to/SRR13340516.0.1x/SRR13340516.0.1x.sorted.merged.MarkDups.BQSR.bam
#/path/to/SRR13340517.0.1x/SRR13340517.0.1x.sorted.merged.MarkDups.BQSR.bam
#/path/to/SRR13340519.0.1x/SRR13340519.0.1x.sorted.merged.MarkDups.BQSR.bam
#/path/to/SRR13340526.0.1x/SRR13340526.0.1x.sorted.merged.MarkDups.BQSR.bam
#/path/to/SRR13340527.0.1x/SRR13340527.0.1x.sorted.merged.MarkDups.BQSR.bam
#/path/to/SRR13340530.0.1x/SRR13340530.0.1x.sorted.merged.MarkDups.BQSR.bam
#/path/to/SRR13340531.0.1x/SRR13340531.0.1x.sorted.merged.MarkDups.BQSR.bam
#/path/to/SRR13340532.0.1x/SRR13340532.0.1x.sorted.merged.MarkDups.BQSR.bam

#obs I missed sample SRR13340529 - so starting with 11 samples ended up with 10 downsampled samples.

###Step 3. Impute chunklists:
2. in ChunkList folders (within 0_scripts) impute the chunks

ChunkList.Chr1 
ChunkList.Chr2 
ChunkList.Chr3-10
ChunkList.Chr11-20
ChunkList.Chr21-30
ChunkList.Chr31-37
ChunkList.Chr38

#e.g.:
mkdir /path/to/4_ImputeAndPhase/ChunkList.Chr38/Down_0.1x
mkdir /path/to/4_ImputeAndPhase/ChunkList.Chr38/Down_0.2x
mkdir /path/to/4_ImputeAndPhase/ChunkList.Chr38/Down_0.3x
mkdir /path/to/4_ImputeAndPhase/ChunkList.Chr38/Down_0.4x
mkdir /path/to/4_ImputeAndPhase/ChunkList.Chr38/Down_0.5x
mkdir /path/to/4_ImputeAndPhase/ChunkList.Chr38/Down_0.6x
mkdir /path/to/4_ImputeAndPhase/ChunkList.Chr38/Down_0.7x
mkdir /path/to/4_ImputeAndPhase/ChunkList.Chr38/Down_0.8x
mkdir /path/to/4_ImputeAndPhase/ChunkList.Chr38/Down_0.9x
mkdir /path/to/4_ImputeAndPhase/ChunkList.Chr38/Down_1x

#Then run 03_ImputeAndPhaseChrsByChunks.sh


###Step 4. Ligate

#chr1
cd /path/to/4_ImputeAndPhase/ChunkList.Chr1/Down_0.1x/
ls *chr1*.bcf > list.chr1.txt
cd /path/to/4_ImputeAndPhase/ChunkList.Chr1/Down_0.2x/
ls *chr1*.bcf > list.chr1.txt
cd /path/to/4_ImputeAndPhase/ChunkList.Chr1/Down_0.3x/
ls *chr1*.bcf > list.chr1.txt
cd /path/to/4_ImputeAndPhase/ChunkList.Chr1/Down_0.4x/
ls *chr1*.bcf > list.chr1.txt
cd /path/to/4_ImputeAndPhase/ChunkList.Chr1/Down_0.5x/
ls *chr1*.bcf > list.chr1.txt
cd /path/to/4_ImputeAndPhase/ChunkList.Chr1/Down_0.6x/
ls *chr1*.bcf > list.chr1.txt
cd /path/to/4_ImputeAndPhase/ChunkList.Chr1/Down_0.7x/
ls *chr1*.bcf > list.chr1.txt
cd /path/to/4_ImputeAndPhase/ChunkList.Chr1/Down_0.8x/
ls *chr1*.bcf > list.chr1.txt
cd /path/to/4_ImputeAndPhase/ChunkList.Chr1/Down_0.9x/
ls *chr1*.bcf > list.chr1.txt
cd /path/to/4_ImputeAndPhase/ChunkList.Chr1/Down_1x/
ls *chr1*.bcf > list.chr1.txt
#etc for all ChunkLists

#in scripts folder:
#create: downsampling.txt
Down_0.1x
Down_0.2x
Down_0.3x
Down_0.4x
Down_0.5x
Down_0.6x
Down_0.7x
Down_0.8x
Down_0.9x
Down_1x

#run script like this: 
while read -r a;do
sbatch LigatePerChromsomeChunkList.DownS.Chr38.sh $a
done < downsampling.txt 


#i.e. slightly updated LigatePerChromsomeChunkList.DownS.Chr1.sh #do one for each:
#!/bin/bash -l
#SBATCH -A $proj_id
#SBATCH -p core -n 1
#SBATCH -J LigatePerChromsomeChunkList.DownS.Chr38
#SBATCH -t 10:00:00

module load bioinfo-tools
module load bcftools
module load samtools

down=$1

Reference='/path/to/UU_Cfam_GSD_1.0_ROSY.fa'
PhaseChunk='/path/to/4_ImputeAndPhase/ChunkList.Chr38'/$down
LigateCHRs='/path/to/5_LigatedImputedCHRs/'
ScriptsFolder='/path/to//path/to/0_Scripts'
Glimpse='/path/to/GLIMPSEv1.1'

cd $PhaseChunk
mkdir $LigateCHRs/$down

# Note! The below script makes an index for the bcf
$Glimpse/GLIMPSE_ligate_static --input list.chr38.txt --output $LigateCHRs/$down/$down.chr38.merged.bcf

bcftools view $LigateCHRs/$down/$down.chr38.merged.bcf -Oz -o $LigateCHRs/$down/$down.chr38.merged.vcf.gz
bcftools index $LigateCHRs/$down/$down.chr38.merged.vcf.gz




##Step 5: QC-filter

#run script like this, one for each downsample-group: 
while read -r a;do
sbatch QCfilter_imputed_chr_Down_0.1x.sh $a
done < CHR1-38.txt


#!/bin/bash -l
#SBATCH -A $proj_id
#SBATCH -p core -n 1
#SBATCH -J QCfilter_imputed_chr
#SBATCH -t 10:00:00
#SBATCH --mail-user $user_email
#SBATCH --mail-type=ALL

module load bioinfo-tools
module load bcftools
module load samtools
module load plink

chrN=$1

LIGATED_IMPUTED_PATH='/path/to/5_LigatedImputedCHRs/Down_0.1x' 
FINAL_IMPUTED_PATH='/path/to/6_final_imputed_dataset' 

cd $LIGATED_IMPUTED_PATH

#qc filter on 
bcftools filter -e 'INFO/INFO<0.8' $LIGATED_IMPUTED_PATH'/Down_0.1x.'$chrN'.merged.vcf.gz' -Oz -o $LIGATED_IMPUTED_PATH'/Down_0.1x.'$chrN'.merged.qc.vcf.gz'

#create SNP ID in gencove:
bcftools annotate --set-id '%CHROM:%POS' $LIGATED_IMPUTED_PATH'/Down_0.1x.'$chrN'.merged.qc.vcf.gz' -Oz -o $LIGATED_IMPUTED_PATH'/Down_0.1x.'$chrN'.merged.qc_modi.vcf.gz'
bcftools index $LIGATED_IMPUTED_PATH'/Down_0.1x.'$chrN'.merged.qc_modi.vcf.gz'


#first create plink files from vcf for each chr #done by sbatch script:
plink --dog --vcf $LIGATED_IMPUTED_PATH'/Down_0.1x.'$chrN'.merged.qc_modi.vcf.gz' --keep-allele-order --double-id --out $FINAL_IMPUTED_PATH'/Down_0.1x.'$chrN --make-bed #--double-id  causes both family and within-family IDs to be set to the sample ID



########################################

FINAL_GENCOVE_PATH=/proj/snic2022-6-229/nobackup/private/Jennifer/DarwinsArkMapping/6_final_imputed_dataset

#make a list to merge: mylist_Down_0.1x_chr1-38.txt
Down_0.1x.chr1
Down_0.1x.chr2
Down_0.1x.chr3
Down_0.1x.chr4
Down_0.1x.chr5
Down_0.1x.chr6




#merge plinkfiles from all chr using the list:

merge_plinkfiles.sh


#!/bin/bash -l
#SBATCH -A $proj_id
#SBATCH -p core -n 1
#SBATCH -J merge_plinkfiles
#SBATCH -t 10:00:00
#SBATCH --mail-user $user_email
#SBATCH --mail-type=ALL

module load bioinfo-tools
module load plink

cd /proj/snic2022-6-229/nobackup/private/Jennifer/DarwinsArkMapping/6_final_imputed_dataset/

plink --allow-no-sex  --dog --geno 0.05 --keep-allele-order --bfile Down_0.1x.chr1 --merge-list mylist_Down_0.1x_chr1-38.txt --make-bed --out Down_0.1x_IMP_ALLCHR
plink --allow-no-sex  --dog --geno 0.05 --keep-allele-order --bfile Down_0.1x.chr1 --merge-list mylist_Down_0.2x_chr1-38.txt --make-bed --out Down_0.2x_IMP_ALLCHR
plink --allow-no-sex  --dog --geno 0.05 --keep-allele-order --bfile Down_0.1x.chr1 --merge-list mylist_Down_0.3x_chr1-38.txt --make-bed --out Down_0.3x_IMP_ALLCHR
plink --allow-no-sex  --dog --geno 0.05 --keep-allele-order --bfile Down_0.1x.chr1 --merge-list mylist_Down_0.4x_chr1-38.txt --make-bed --out Down_0.4x_IMP_ALLCHR
plink --allow-no-sex  --dog --geno 0.05 --keep-allele-order --bfile Down_0.1x.chr1 --merge-list mylist_Down_0.5x_chr1-38.txt --make-bed --out Down_0.5x_IMP_ALLCHR
plink --allow-no-sex  --dog --geno 0.05 --keep-allele-order --bfile Down_0.1x.chr1 --merge-list mylist_Down_0.6x_chr1-38.txt --make-bed --out Down_0.6x_IMP_ALLCHR
plink --allow-no-sex  --dog --geno 0.05 --keep-allele-order --bfile Down_0.1x.chr1 --merge-list mylist_Down_0.7x_chr1-38.txt --make-bed --out Down_0.7x_IMP_ALLCHR
plink --allow-no-sex  --dog --geno 0.05 --keep-allele-order --bfile Down_0.1x.chr1 --merge-list mylist_Down_0.8x_chr1-38.txt --make-bed --out Down_0.8x_IMP_ALLCHR
plink --allow-no-sex  --dog --geno 0.05 --keep-allele-order --bfile Down_0.1x.chr1 --merge-list mylist_Down_0.9x_chr1-38.txt --make-bed --out Down_0.9x_IMP_ALLCHR
plink --allow-no-sex  --dog --geno 0.05 --keep-allele-order --bfile Down_0.1x.chr1 --merge-list mylist_Down_1x_chr1-38.txt --make-bed --out Down_1x_IMP_ALLCHR


plink2 --dog --bfile Down_0.1x_IMP_ALLCHR --maf 0.05 --make-bed --out Down_0.1x_IMP_ALLCHR.maf0.05
#6012422 variants remaining after main filters.
plink2 --dog --bfile Down_0.2x_IMP_ALLCHR --maf 0.05 --make-bed --out Down_0.2x_IMP_ALLCHR.maf0.05
#7466301 variants remaining after main filters.
plink2 --dog --bfile Down_0.3x_IMP_ALLCHR --maf 0.05 --make-bed --out Down_0.3x_IMP_ALLCHR.maf0.05
#7730501 variants remaining after main filters.
plink2 --dog --bfile Down_0.4x_IMP_ALLCHR --maf 0.05 --make-bed --out Down_0.4x_IMP_ALLCHR.maf0.05
#7817142 variants remaining after main filters.
plink2 --dog --bfile Down_0.5x_IMP_ALLCHR --maf 0.05 --make-bed --out Down_0.5x_IMP_ALLCHR.maf0.05
#7857791 variants remaining after main filters.
plink2 --dog --bfile Down_0.6x_IMP_ALLCHR --maf 0.05 --make-bed --out Down_0.6x_IMP_ALLCHR.maf0.05
#7880800 variants remaining after main filters.
plink2 --dog --bfile Down_0.7x_IMP_ALLCHR --maf 0.05 --make-bed --out Down_0.7x_IMP_ALLCHR.maf0.05
#7895854 variants remaining after main filters.
plink2 --dog --bfile Down_0.8x_IMP_ALLCHR --maf 0.05 --make-bed --out Down_0.8x_IMP_ALLCHR.maf0.05
plink2 --dog --bfile Down_0.9x_IMP_ALLCHR --maf 0.05 --make-bed --out Down_0.9x_IMP_ALLCHR.maf0.05
plink2 --dog --bfile Down_1x_IMP_ALLCHR --maf 0.05 --make-bed --out Down_1x_IMP_ALLCHR.maf0.05


#check number of variants after maf filter for each set: 0.1x etc
#make list of dogs
plink2 --dog --bfile Down_ALLx_IMP_ALLCHR --maf 0.01 --keep --make-bed --out Down_ALLx_IMP_ALLCHR.maf0.01

#obs I first ran for the maf filted data:
plink2 --dog --bfile Down_0.1x_IMP_ALLCHR.maf0.01 --out sample_counts_Down_0.1x_IMP --sample-counts
#etc

#then for non-filtered data:
plink2 --dog --bfile Down_0.1x_IMP_ALLCHR --out allsample_counts_Down_0.1x_IMP --sample-counts
#etc

####check depth - calculate with

/path/to/SRR13340532.0.9x/SRR13340532.0.9x.sorted.merged.MarkDups.BQSR.bam.knownsites.depth.txt

#make script: statMove.sh
chmod +x statMove.sh

#!/bin/bash
# To run
#               ./statMove.sh Darwin_DOWNstats.movelist      #
#

while read line1
do

######## Variables to change #########
dog_ID=$line1

# Path to folder to store stat files
outputFolder='/path/to/mapping_stats_DOWN'
# Path to output folder
baseFolder='/path/to'

echo "#!/bin/bash -l" > statMove.sc;
echo "#SBATCH -A naiss2024-5-494" >> statMove.sc;
echo "#SBATCH -p core -n 1" >> statMove.sc;
echo "#SBATCH -J $dog_ID-bamMove" >> statMove.sc;
echo "#SBATCH -t 10:00" >> statMove.sc;

echo "dog_ID=$dog_ID" >> statMove.sc;
echo "outputFolder=$outputFolder" >> statMove.sc;
echo "baseFolder=$baseFolder" >> statMove.sc;

echo "cd $baseFolder/$dog_ID" >> statMove.sc;
echo "cp $dog_ID.sorted.merged.MarkDups.BQSR.bam.knownsites.depth.txt $outputFolder" >> statMove.sc;

sbatch statMove.sc;

done<Darwin_DOWNstats.movelist 

###########################

./extract_depth_autosomes_and_x.sh

#!/bin/bash

output_file="combined_knownsites.depth.txt" #found in data/03_qc

# Clear the output file if it exists
> $output_file

# Use a pattern that matches SRR followed by numbers and letters at the start of filenames
for file in SRR*.txt; do
  # Extract sample ID from filename (entire name before .txt)
  sample_id=$(basename "$file" .txt)
  
  # Extract 3rd value from 2nd and 3rd rows
  values=$(awk 'NR==2 || NR==3 {print $3}' "$file")
  
  # Combine sample ID with extracted values and append to output file
  echo -e "${sample_id}\t$(echo $values | tr '\n' '\t')" >> $output_file
done




#######################################
#next script: 3_genotype_concordance.sh
######################################

