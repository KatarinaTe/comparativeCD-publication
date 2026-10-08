#!/bin/bash
# To run
#               ./02_GenotypeLiklihood.sh CHR1-38.txt      #
# Note! Extract variable SNP sites per chr from phased Dog10K reference panel

while read line1
do

## Variables to change ##

CHR=$line1

ScriptsFolder='/path/to/0_Scripts'
VCFandTSVpath='/path/to/1_phased-imputation-panel'
REFGEN='/path/to/UU_Cfam_GSD_1.0_ROSY.fa'
OUT='/path/to/3_GenotypeLikelihood'


echo "#!/bin/bash -l" > 02_GenotypeLiklihood.sc;
echo "#SBATCH -A naiss2023-5-396" >> 02_GenotypeLiklihood.sc;
echo "#SBATCH -p core -n 1" >> 02_GenotypeLiklihood.sc;
echo "#SBATCH -J 02_GenotypeLiklihood.$CHR" >> 02_GenotypeLiklihood.sc;
echo "#SBATCH -t 96:00:00" >> 02_GenotypeLiklihood.sc;

## load some modules ## 
echo "module load bioinfo-tools" >> 02_GenotypeLiklihood.sc;
echo "module load bcftools" >> 02_GenotypeLiklihood.sc;
echo "module load samtools" >> 02_GenotypeLiklihood.sc;

echo "CHR=$CHR" >> 02_GenotypeLiklihood.sc;
echo "ScriptsFolder=$ScriptsFolder" >> 02_GenotypeLiklihood.sc;
echo "VCFandTSVpath=$VCFandTSVpath" >> 02_GenotypeLiklihood.sc;
echo "REFGEN=$REFGEN" >> 02_GenotypeLiklihood.sc;
echo "OUT=$OUT" >> 02_GenotypeLiklihood.sc;


echo "cd ${OUT}" >> 02_GenotypeLiklihood.sc;

 
echo "bcftools mpileup -f ${REFGEN} -I -E -a 'FORMAT/DP' -T ${VCFandTSVpath}/${CHR}.Dog10K.Phased.snp.sites.vcf.gz -r ${CHR} --bam-list ${ScriptsFolder}/All.DA_LowPass.Samples_bam.txt -Ou | bcftools call -Aim -C alleles -T ${VCFandTSVpath}/${CHR}.Dog10K.Phased.snp.sites.tsv.gz -Oz -o ${OUT}/${CHR}.GL.DA_LowPass.vcf.gz" >> 02_GenotypeLiklihood.sc;
echo "bcftools index -f ${OUT}/${CHR}.GL.DA_LowPass.vcf.gz" >> 02_GenotypeLiklihood.sc;
sbatch 02_GenotypeLiklihood.sc;

done<CHR1-38.txt