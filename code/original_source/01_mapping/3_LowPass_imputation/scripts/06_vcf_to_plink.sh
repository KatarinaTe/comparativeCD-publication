#!/bin/bash -l
#SBATCH -A $proj_id
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 24:00:00
#SBATCH -J vcf_to_plink
#SBATCH --mail-user $user_email
#SBATCH --mail-type=ALL

#run script like this: 
#while read -r a;do
#sbatch sbatch_vcf_to_plink.sh $a
#done < CHR1-38.txt 

chrN=$1

module load bioinfo-tools plink 

GENCOVE_PATH=/path/to/5_LigatedImputedCHRs/ #gencove imputed are here

cd $GENCOVE_PATH

#first create plink files from vcf for each chr
plink --dog --vcf 'DA.'$chrN'.merged.qc_modi.vcf.gz' --keep-allele-order --double-id --out 'DA_IMP_GENCOVE_'$chrN --make-bed #--double-id  causes both family and within-family IDs to be set to the sample ID
