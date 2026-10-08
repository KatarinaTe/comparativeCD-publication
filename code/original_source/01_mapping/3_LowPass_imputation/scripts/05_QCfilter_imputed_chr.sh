#!/bin/bash
# To run
#               ./05_QCfilter_imputed_chr.sh CHR1-38.txt      #


while read line1
do

## Variables to change ##

chrN=$line1


echo "#!/bin/bash -l" > 05_QCfilter_imputed_chr.sc;
echo "#SBATCH -A $proj_id" >> 05_QCfilter_imputed_chr.sc;
echo "#SBATCH -p core -n 1" >> 05_QCfilter_imputed_chr.sc;
echo "#SBATCH -J 05_QCfilter$chrN" >> 05_QCfilter_imputed_chr.sc;
echo "#SBATCH -t 10:00:00" >> 05_QCfilter_imputed_chr.sc;

## load some modules ## 
echo "module load bioinfo-tools" >> 05_QCfilter_imputed_chr.sc;
echo "module load bcftools" >> 05_QCfilter_imputed_chr.sc;
echo "module load samtools" >> 05_QCfilter_imputed_chr.sc;

echo "cd /path/to/5_LigatedImputedCHRs" >> 05_QCfilter_imputed_chr.sc;

echo "chrN=$chrN" >> 05_QCfilter_imputed_chr.sc;
echo "bcftools filter -e 'INFO/INFO<0.8' DA.$chrN.merged.vcf.gz -Oz -o DA.$chrN.merged.qc.vcf.gz" >> 05_QCfilter_imputed_chr.sc;
echo "bcftools annotate --set-id '%CHROM:%POS' DA.$chrN.merged.qc.vcf.gz -Oz -o DA.$chrN.merged.qc_modi.vcf.gz" >> 05_QCfilter_imputed_chr.sc;

sbatch 05_QCfilter_imputed_chr.sc;

done<CHR1-38.txt
