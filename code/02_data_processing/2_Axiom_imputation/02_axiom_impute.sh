################################################
################################################
#Imputation to dog10K!
################################################
################################################

DA_AFFY #lifted to canfam4: consists of 411 individuals (specified in ind_to_keep_round3.txt) that we use for the Darwins Ark dataset to be merged with the rest (gencove).

##############################################################
##run all chr with conform.gt and beagle: axiom_impute.sh:
##############################################################

#!/bin/bash -l
#SBATCH -A $proj_id
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 2-00:00:00
#SBATCH -J imputation_axiom
#SBATCH --mail-user $user_email
#SBATCH --mail-type=ALL

module load bioinfo-tools tabix bcftools plink2

#run script like this: 
#while read -r a b c;do
#sbatch axiom_impute.sh $a $b $c
#done < cf4_chr_start_stop.bed #in ./data folder

chrN=$1
startN=$2
stopN=$3

RUN_PATH=/path/to/my_folder/
REF_PATH=/path/to/ref_folder/
cd $RUN_PATH


#create vcf file for each chr
plink2 --bfile DA_AFFY --dog --chr $chrN --recode vcf id-paste=iid  --fa $REF_PATH'UU_Cfam_GSD_1.0_ROSY.fa' --ref-from-fa --out DA_AFFY_$chrN

#make vcf.gz and index
bcftools view DA_AFFY_$chrN.vcf  -Oz -o DA_AFFY_$chrN.vcf.gz 
bcftools index DA_AFFY_$chrN.vcf.gz


#rename the chr ID using the existing file chromosomes.txt
bcftools annotate --rename-chrs chromosomes.txt 'DA_AFFY_'$chrN'.vcf' -Oz -o 'DA_AFFY_'$chrN'_rename.vcf.gz'
bcftools index 'DA_AFFY_'$chrN'_rename.vcf.gz'

#make one ref per chr: DO THIS ONCE
#bcftools view /path/to/AutoAndXPAR.Dog10K.Phased.bcf -r $chrN  -Oz -o $REF_PATH'modi_AutoAndXPAR.Dog10K.Phased_'$chrN'.vcf.gz'
#bcftools index $REF_PATH/'modi_AutoAndXPAR.Dog10K.Phased_'$chrN'.vcf.gz'

### run conform gt ###strict=TRUE
java -jar ../conform-gt.24May16.cee.jar ref=$REF_PATH'modi_AutoAndXPAR.Dog10K.Phased_'$chrN'.vcf.gz' gt='DA_AFFY_'$chrN'_rename.vcf.gz' chrom=$chrN match=POS out='DA_AFFY_conf_'$chrN strict=TRUE
bcftools index 'DA_AFFY_conf_'$chrN'.vcf.gz'

### run beagle ### iterations=40 burnin=30
java -Xmx50g -jar ../beagle.22Jul22.46e.jar gt='DA_AFFY_conf_'$chrN'.vcf.gz' chrom=$chrN":"$startN"-"$stopN map=/data/mapfiles/'canFam4.cM.'$chrN'.map' out='DA_AFFYimp_'$chrN  ref=$REF_PATH'modi_AutoAndXPAR.Dog10K.Phased_'$chrN'.vcf.gz' iterations=40 burnin=30
bcftools index 'DA_AFFYimp_'$chrN'.vcf.gz'


## I had settings: java –Xmx4g: But it actually ran: Command line: java -Xmx3641m -jar beagle.22Jul22.46e.jar
#nthreads= If the nthreads parameter is not specified, the nthreads parameter will be set equal to the number of CPU cores on the host machine.
#it ran: nthreads=20
#but need to increase: java -Xmx50g -jar beagle.22Jul22.46e.jar (needed for chr1)

#quality filter the output from beagle
bcftools filter -e 'INFO/DR2<0.8' 'DA_AFFYimp_'$chrN'.vcf.gz' -Oz -o 'DA_AFFYimp_'$chrN'.qc.vcf.gz'
bcftools index 'DA_AFFYimp_'$chrN'.qc.vcf.gz'

#give names to the SNPs
bcftools annotate --set-id '%CHROM:%POS' 'DA_AFFYimp_'$chrN'.qc.vcf.gz' -Oz -o 'DA_AFFYimp_'$chrN'.qc.modi.vcf.gz'
bcftools index 'DA_AFFYimp_'$chrN'.qc.modi.vcf.gz'

#check number of variants and individuals:
echo 'bcftools query -l 'DA_AFFYimp_'$chrN'.qc.modi.vcf.gz' | wc -l '
bcftools query -l 'DA_AFFYimp_'$chrN'.qc.modi.vcf.gz' | wc -l 

echo 'bcftools view -H 'DA_AFFYimp_'$chrN'.vcf.gz' | wc -l'
bcftools view -H 'DA_AFFYimp_'$chrN'.vcf.gz' | wc -l

echo 'bcftools view -H 'DA_AFFYimp_'$chrN'.qc.modi.vcf.gz'| wc -l'
bcftools view -H 'DA_AFFYimp_'$chrN'.qc.modi.vcf.gz'| wc -l

