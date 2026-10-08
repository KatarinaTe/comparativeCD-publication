#!/bin/bash -l
#SBATCH -A $proj_id
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 12:00:00
#SBATCH -J merge_chr_plink
#SBATCH --mail-user $user_email
#SBATCH --mail-type=ALL

module load bioinfo-tools plink 

GENCOVE_PATH=/path/to/5_LigatedImputedCHRs/ #gencove imputed are here
FINAL_GENCOVE_PATH=/path/to/6_Final_imputed_dataset

#make this list according to your samples names: mylist_DA_IMP_GENCOVE_chr1-38.txt
#DA_IMP_GENCOVE_chr1
#DA_IMP_GENCOVE_chr2
#DA_IMP_GENCOVE_chr3 etc.

cd $GENCOVE_PATH

plink --allow-no-sex  --dog --geno 0.05 --keep-allele-order --bfile DA_IMP_GENCOVE_chr1 --merge-list mylist_DA_IMP_GENCOVE_chr1-38.txt --make-bed --out DA_IMP_GENCOVE_ALLCHR