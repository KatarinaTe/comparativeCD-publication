####################################################################
#AXIOM create plink files out of vcf and then merge chr using plink:
####################################################################

###############################################################
#1. create plink files from vcf for each chr using batch script
###############################################################

#run script like this: 
while read -r a;do
sbatch sbatch_vcf_to_plink_axiomONLY.sh $a
done < chr.txt # i.e. a list of chromosomes:
#chr1
#chr2
#etc

###start script:
#!/bin/bash -l
#SBATCH -A $proj_id
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 24:00:00
#SBATCH -J imp_axiomvcf_to_plink
#SBATCH --mail-type=ALL

chrN=$1

module load bioinfo-tools plink

RUN_PATH=/path/to/run_folder
FINAL_IMPFILE_PATH=/path/to/final_imp_folder

cd $RUN_PATH

#first create plink files from vcf for each chr
plink --dog --vcf 'DA_AFFYimp_'$chrN'.qc.modi.vcf.gz' --keep-allele-order --double-id --out $FINAL_IMPFILE_PATH'/DA_AFFYimp_'$chrN --make-bed 
#--double-id  causes both family and within-family IDs to be set to the sample ID
####################################################
#2. Merge chromosomes using plink:
#####################################################

#make a list to merge: mylist_DA_AFFYimp_chr1-38.txt #
#DA_AFFYimp_chr1.qc.modi.vcf.gz
#DA_AFFYimp_chr2.qc.modi.vcf.gz
#DA_AFFYimp_chr3.qc.modi.vcf.gz
#etc

#merge plinkfiles from all chr using the list:
plink --allow-no-sex  --dog --geno 0.05 --keep-allele-order --bfile DA_AFFYimp_chr1 --merge-list mylist_DA_AFFYimp_chr1-38.txt --make-bed --out DA_AFFYimp_ALLCHR
#14042437 variants and 411 dogs pass filters and QC.