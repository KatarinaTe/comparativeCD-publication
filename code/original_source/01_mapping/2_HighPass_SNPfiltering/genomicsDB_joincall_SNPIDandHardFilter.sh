#! /bin/bash
#SBATCH -A $proj_id
#SBATCH -p core
#SBATCH -n 10
#SBATCH -t 4-00:00:00
#SBATCH -J Join_calling
#SBATCH --mail-user $user_email
#SBATCH --mail-type=ALL


##load the modules
module load bioinfo-tools
module load bwa/0.7.17
module load samtools
module load bcftools
module load tabix
#module load picard
module load GATK/4.2.0.0
module load picard/3.1.1 java/OpenJDK_17+35 ##new edit 2024-09-20!!!


## set the path for ref and others
#reads_path=/path/to/samples
assembly_ref=/path/to/UU_Cfam_GSD_1.0_ROSY.fa
assembly_base=$(basename $assembly_ref)
#sample=($(awk '{print $1}' $1))
#sample_id=${sample[0]}
#runs=($(awk '{print $2}' $1))
known_variants=/path/to/UU_Cfam_GSD_1.0.BQSR.DB.bed.gz
chunks_path=/path/to/chunks_10/
#bwa2_path=/path/to/bwa-mem2-2.1_x64-linux/

##uses cohort.sample_map in ./data folder put this in the run folder: /path/to/3_Joint_calling


chunk_process_s1(){
	i=$1
	chunk_i=$(basename $i ".bed")

	gatk --java-options "-Xmx4G" GenomicsDBImport \
		--genomicsdb-workspace-path 'my_database.'$chunk_i \
		--sample-name-map cohort.sample_map \
		--batch-size 50 \
		--genomicsdb-shared-posixfs-optimizations true \
		-L $i

	gatk --java-options "-Xmx4G" GenotypeGVCFs \
	    -R $assembly_ref \
	    -V gendb://'my_database.'$chunk_i \
	    -L $i \
	    -O $chunk_i'.raw_output.vcf.gz'
}

chunk_process_s2(){
	chunk_i='chunk_'$1
	tabix --list-chroms $chunk_i'.raw_output.vcf.gz' > $chunk_i'.list'
	zcat $chunk_i'.raw_output.vcf.gz' |grep "^#" > $chunk_i'.header'
	while IFS= read -r line;do
		tabix $chunk_i'.raw_output.vcf.gz' $line |cat $chunk_i'.header' - |bgzip > $line'.raw_output.vcf.gz'
		tabix -p vcf $line'.raw_output.vcf.gz'
	done < $chunk_i'.list'
}

for j in $chunks_path/*.bed;do
	(
	chunk_process_s1 $j
	)&
done
wait

for j in {1..10};do
	(
	chunk_process_s2 $j
	)&
done
wait


gatk --java-options "-Xmx6G" MergeVcfs \
-I chr1.raw_output.vcf.gz \
-I chr2.raw_output.vcf.gz \
-I chr3.raw_output.vcf.gz \
-I chr4.raw_output.vcf.gz \
-I chr5.raw_output.vcf.gz \
-I chr6.raw_output.vcf.gz \
-I chr7.raw_output.vcf.gz \
-I chr8.raw_output.vcf.gz \
-I chr9.raw_output.vcf.gz \
-I chr10.raw_output.vcf.gz \
-I chr11.raw_output.vcf.gz \
-I chr12.raw_output.vcf.gz \
-I chr13.raw_output.vcf.gz \
-I chr14.raw_output.vcf.gz \
-I chr15.raw_output.vcf.gz \
-I chr16.raw_output.vcf.gz \
-I chr17.raw_output.vcf.gz \
-I chr18.raw_output.vcf.gz \
-I chr19.raw_output.vcf.gz \
-I chr20.raw_output.vcf.gz \
-I chr21.raw_output.vcf.gz \
-I chr22.raw_output.vcf.gz \
-I chr23.raw_output.vcf.gz \
-I chr24.raw_output.vcf.gz \
-I chr25.raw_output.vcf.gz \
-I chr26.raw_output.vcf.gz \
-I chr27.raw_output.vcf.gz \
-I chr28.raw_output.vcf.gz \
-I chr29.raw_output.vcf.gz \
-I chr30.raw_output.vcf.gz \
-I chr31.raw_output.vcf.gz \
-I chr32.raw_output.vcf.gz \
-I chr33.raw_output.vcf.gz \
-I chr34.raw_output.vcf.gz \
-I chr35.raw_output.vcf.gz \
-I chr36.raw_output.vcf.gz \
-I chr37.raw_output.vcf.gz \
-I chr38.raw_output.vcf.gz \
-I chrX.raw_output.vcf.gz \
-I chrY_NC_051844.1.raw_output.vcf.gz \
-I chrM.raw_output.vcf.gz \
-I chrY_unplaced_NW_024010443.1.raw_output.vcf.gz \
-I chrY_unplaced_NW_024010444.1.raw_output.vcf.gz \
-I chunk_11.raw_output.vcf.gz \
-O raw_output2.vcf.gz &&

gatk --java-options "-Xmx5g -Xms5g" VariantRecalibrator \
-R $assembly_ref \
-V raw_output.vcf.gz \
--resource:array,known=false,training=true,truth=true,prior=12.0 /proj/snic2021-6-208/jennifer/F.G.vanSteenbeek/3_Joint_calling/SRZ189891_722g.simp.header.CanineHDandAxiom_K9_HD.GSD_1.0.vcf.gz \
-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP \
-mode SNP \
-O snp.output.recal \
--tranches-file snp.output.tranches \
--rscript-file snp.output.plots.R \
--trust-all-polymorphic true

#gatk --java-options "-Xmx5g -Xms5g" ApplyVQSR \
#-V raw_output.vcf.gz \
#--recal-file snp.output.recal \
#--tranches-file snp.output.tranches \
#--truth-sensitivity-filter-level 99.0 \
#--create-output-variant-index true \
#-mode SNP \
#-O snp.recal99.vcf.gz

####SNPIDandFilterSteps:
# STEP 1 ## Annotate SNPs
bcftools annotate --set-id '%CHROM\:%POS' raw_output2.vcf.gz Oz -o raw_output.ann.id.vcf.gz
tabix raw_output.ann.id.vcf.gz

## STEP 4 ## Hard filters
# Run HF SNP on All - as per GATK 
gatk SelectVariants  -V raw_output.ann.id.vcf.gz  -select-type SNP -O SNP.raw_output.ann.id.vcf.gz

gatk VariantFiltration -V SNP.raw_output.ann.id.vcf.gz -filter "QD < 2.0" --filter-name "QD2" -filter "QUAL < 30.0" --filter-name "QUAL30" -filter "SOR > 3.0" --filter-name "SOR3" -filter "FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" -O SNP.HF.ann.id.vcf.gz
tabix SNP.HF.ann.id.vcf.gz

#zcat SNP.HF.ann.id.vcf.gz| grep -v "#" | less -S
