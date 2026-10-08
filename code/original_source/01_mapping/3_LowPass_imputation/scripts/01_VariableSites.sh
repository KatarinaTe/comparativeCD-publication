#!/bin/bash
# To run
#               ./01_VariableSites.sh CHR1-38.txt      #
# Note! Extract variable SNP sites per chr from phased Dog10K reference panel

while read line1
do

## Variables to change ##

CHR=$line1

ScriptsFolder='/path/to/0_Scripts'
BCF='/path/to/1_phased-imputation-panel/AutoAndXPAR.Dog10K.Phased.bcf'
Reference='/path/to/UU_Cfam_GSD_1.0_ROSY.fa'
OUT='/path/to/1_phased-imputation-panel'


echo "#!/bin/bash -l" > 01_VariableSites.sc;
echo "#SBATCH -A $proj_id" > 01_VariableSites.sc;
echo "#SBATCH -p core -n 1" > 01_VariableSites.sc;
echo "#SBATCH -J VariableSites_Dog10KReferencePanel.$CHR" > 01_VariableSites.sc;
echo "#SBATCH -t 24:00:00" > 01_VariableSites.sc;
echo "#SBATCH -o /path/to/Z_ErrorAndOut/VariableSites_Dog10KReferencePanel.$CHR.output" > 01_VariableSites.sc;
echo "#SBATCH -e /path/to/Z_ErrorAndOut/VariableSites_Dog10KReferencePanel.$CHR.error" > 01_VariableSites.sc;

## load some modules ## 
echo "module load bioinfo-tools" > 01_VariableSites.sc;
echo "module load bcftools" > 01_VariableSites.sc;
echo "module load samtools" > 01_VariableSites.sc;

echo "CHR=$CHR" > 01_VariableSites.sc;
echo "ScriptsFolder=$ScriptsFolder" > 01_VariableSites.sc;
echo "BCF=$BCF" > 01_VariableSites.sc;
echo "Reference=$Reference" > 01_VariableSites.sc;
echo "OUT=$OUT" > 01_VariableSites.sc;


echo "cd $OUT" > 01_VariableSites.sc;

# We a VCF/BCF file containing only sites to tell BCFtools at which positions to make a genotype call. 
# Since BCFtools does not compute correctly genotype likelihood for indels, here we only focus on SNPs (however, GLIMPSE can impute any type of variants as soon it is bi-allelic and has GLs being properly defined). 
# To perform the extraction from one chromosome of a reference panel, run first BCFtools as follows:

# -G -m 2 -M 2 -v snps, -G group for HWE filter, -m multiallic caller for rare variant calling, -M keep masked, -v snps keep variant sites only

echo "bcftools view -G -m 2 -M 2 -v snps -r $CHR $BCF -Oz -o $CHR.Dog10K.Phased.snp.sites.vcf.gz" > 01_VariableSites.sc;
echo "bcftools index -f $CHR.Dog10K.Phased.snp.sites.vcf.gz" > 01_VariableSites.sc;

# Then, convert the output file $CHR.Dog10K.Phased.snp.sites.vcf.gz into TSV format and index the file using tabix (requires htslib in PATH), using this command:

echo "bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' $CHR.Dog10K.Phased.snp.sites.vcf.gz | bgzip -c > $CHR.Dog10K.Phased.snp.sites.tsv.gz" > 01_VariableSites.sc;
echo "tabix -s1 -b2 -e2 $CHR.Dog10K.Phased.snp.sites.tsv.gz" > 01_VariableSites.sc;

sbatch 01_VariableSites.sc;

done<CHR1-38.txt