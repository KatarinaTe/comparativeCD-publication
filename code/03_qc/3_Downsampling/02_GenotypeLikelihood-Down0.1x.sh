#!/bin/bash
# To run
#               ./02_GenotypeLikelihood-Down0.1x.sh CHR1-38.txt      #
# Note! Extract variable SNP sites per chr from phased Dog10K reference panel

while read line1
do

## Variables to change ##

CHR=$line1

ScriptsFolder='/proj/snic2022-6-229/nobackup/private/Jennifer/DarwinsArkMapping/0_Scripts'
VCFandTSVpath='/proj/snic2022-6-229/nobackup/private/Jennifer/DarwinsArkMapping/1_phased-imputation-panel'

REFGEN='/crex/proj/uppstore2017228/KLT.03.CGEN/Mischka/Other_dog_genome_assembly/GSD1.0_RosY/UU_Cfam_GSD_1.0_ROSY.fa'
OUT='/proj/snic2022-6-229/nobackup/private/Jennifer/DarwinsArkMapping/3_GenotypeLikelihood/3_GenotypeLikelihood_0.1x'


echo "#!/bin/bash -l" > 02_GenotypeLikelihood-Down0.1x.sc;
echo "#SBATCH -A naiss2024-5-494" >> 02_GenotypeLikelihood-Down0.1x.sc;
echo "#SBATCH -p core -n 1" >> 02_GenotypeLikelihood-Down0.1x.sc;
echo "#SBATCH -J 02_GenotypeLiklihood.$CHR" >> 02_GenotypeLikelihood-Down0.1x.sc;
echo "#SBATCH -t 96:00:00" >> 02_GenotypeLikelihood-Down0.1x.sc;

## load some modules ## 
echo "module load bioinfo-tools" >> 02_GenotypeLikelihood-Down0.1x.sc;
echo "module load bcftools" >> 02_GenotypeLikelihood-Down0.1x.sc;
echo "module load samtools" >> 02_GenotypeLikelihood-Down0.1x.sc;

echo "CHR=$CHR" >> 02_GenotypeLikelihood-Down0.1x.sc;
echo "ScriptsFolder=$ScriptsFolder" >> 02_GenotypeLikelihood-Down0.1x.sc;
echo "VCFandTSVpath=$VCFandTSVpath" >> 02_GenotypeLikelihood-Down0.1x.sc;
echo "REFGEN=$REFGEN" >> 02_GenotypeLikelihood-Down0.1x.sc;
echo "OUT=$OUT" >> 02_GenotypeLikelihood-Down0.1x.sc;

echo "cd ${OUT}" >> 02_GenotypeLikelihood-Down0.1x.sc;
 
echo "bcftools mpileup -f ${REFGEN} -I -E -a 'FORMAT/DP' -T ${VCFandTSVpath}/${CHR}.Dog10K.Phased.snp.sites.vcf.gz -r ${CHR} --bam-list ${ScriptsFolder}/DownSampledBams.0.1x.txt -Ou | bcftools call -Aim -C alleles -T ${VCFandTSVpath}/${CHR}.Dog10K.Phased.snp.sites.tsv.gz -Oz -o ${OUT}/${CHR}.GL.DownSampled.0.1x.vcf.gz" >> 02_GenotypeLikelihood-Down0.1x.sc;
echo "bcftools index -f ${OUT}/${CHR}.GL.DownSampled.0.1x.vcf.gz" >> 02_GenotypeLikelihood-Down0.1x.sc;
sbatch 02_GenotypeLikelihood-Down0.1x.sc;

done<CHR1-38.txt

