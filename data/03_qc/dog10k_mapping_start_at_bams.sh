#! /bin/bash
#SBATCH -A $proj_id
#SBATCH -p core
#SBATCH -n 10
#SBATCH -t 2-00:00:00
#SBATCH -J Illumina_mapping
##SBATCH --mail-user $user_email
##SBATCH --mail-type=ALL

##load the modules
module load bioinfo-tools
module load bwa/0.7.17
module load samtools
module load GATK/4.1.4.1
module load python/3.7.2
#module load picard
module load picard/3.1.1 java/OpenJDK_17+35 #new edit 2024-09-20!!!

## set the path for ref and others
#reads_path=/path/to/samples #could be used instead of input file:
input_file=$1

sample=($(awk '{print $1}' $input_file))
sample_id=${sample[0]}
#runs=($(awk '{print $2}' $input_file)) #not needed here

#all of these ref genome files are found here: https://github.com/jmkidd/dogmap
assembly_ref=/path/to/UU_Cfam_GSD_1.0_ROSY.fa
assembly_base=$(basename $assembly_ref)
known_variants=/path/to/UU_Cfam_GSD_1.0.BQSR.DB.bed.gz
depth_sites=/path/to/SRZ189891_722g.simp.header.CanineHD.names.GSD_1.0.filter.vcf.gz

#addtional paths needed (see https://github.com/Chao912/dog_10k/blob/main/dog_10k_mapping.sh )
chunks_path=/path/to/chunks_10/ #split the genome into 10 chunks  
bwa2_path=/path/to/bwa-mem2-2.1_x64-linux/
gatk_path=/path/to/gatk-4.2.0.0/

#available in the 1_High.LowPass_mapping folder:
run_stats=/path/to/run-stats.py

mkdir -p $sample_id && cd $sample_id

##mapping

#for run_id in ${runs[@]};do
while read -r sample_id run_id bam;do
	echo $sample_id
	echo $run_id
	echo $bam
	#echo $r2_fq
	#$bwa2_path'/bwa-mem2' mem -K 100000000  -t 10 -Y \
	#-R '@RG\tID:'$run_id'\tPL:illumina\tLB:'$sample_id'\t'SM:$sample_id \
	#$assembly_ref $r1_fq $r2_fq |samtools view -hbS > $sample_id'_'$run_id'.bam'
	samtools view -hbS $bam -o $sample_id'_'$run_id'.bam'
	samtools sort -@ 10 -m 1G $sample_id'_'$run_id'.bam' -o $sample_id'_'$run_id'.sorted.bam' && rm $sample_id'_'$run_id'.bam'
	samtools index $sample_id'_'$run_id'.sorted.bam'
done < ../$input_file

sample=($(awk '{print $1}' ../$input_file))
sample_id=${sample[0]}

###merge
ls $sample_id'_'*'sorted.bam' > $sample_id'.bamlist'
if [ $(cat $sample_id'.bamlist'|wc -l) -gt 1 ];then
	samtools merge -@ 10 -b $sample_id'.bamlist' $sample_id'.sorted.merged.bam'
else
	mv $sample_id*'sorted.bam' $sample_id'.sorted.merged.bam'
fi &&
xargs -a $sample_id'.bamlist' rm
rm $sample_id'_'*'sorted.bam.bai'
#
samtools index $sample_id'.sorted.merged.bam'
###extract the unmapped reads
samtools view  -f 4 $sample_id'.sorted.merged.bam' |perl -lane '$line=$_;if($F[0]=~/@/){print}else{print if $F[2] eq "*"}' |samtools view -b > $sample_id'.unmapped.bam'

##split the bam files, and generate the BQSR table for each chunk
chunk_process_s1(){
	i=$1
	chunk_i=$(basename $i ".bed")
	##split the bam
	samtools view -hb -L $i $sample_id'.sorted.merged.bam' > $sample_id'.sorted.merged.'$chunk_i'.bam' &&

	##mark duplicates
	java -Xmx3G -jar $PICARD_ROOT/picard.jar MarkDuplicates \
	METRICS_FILE=$sample_id'.sorted.merged.'$chunk_i'.matrix' \
	INPUT=$sample_id'.sorted.merged.'$chunk_i'.bam' \
	OUTPUT=$sample_id'.sorted.merged.MarkDups.'$chunk_i'.bam'

	samtools index $sample_id'.sorted.merged.MarkDups.'$chunk_i'.bam'

	##BQSR
	# Generate the first pass BQSR table file
	$gatk_path/gatk --java-options "-Xmx4G"  BaseRecalibrator \
	-R $assembly_ref \
	-I $sample_id'.sorted.merged.MarkDups.'$chunk_i'.bam' \
	-L $i \
	--known-sites $known_variants \
	-O $sample_id'.sorted.merged.MarkDups.'$chunk_i'.table'
}

## apply the BQSR for each chunk
chunk_process_s2(){
	i=$1
	chunk_i=$(basename $i ".bed")
	# Apply BQSR
	$gatk_path/gatk --java-options "-Xmx4G"  ApplyBQSR \
	-R $assembly_ref \
	-I $sample_id'.sorted.merged.MarkDups.'$chunk_i'.bam' \
	-L $i \
	-bqsr $sample_id'.BQSR.reports.table' \
	--preserve-qscores-less-than 6 \
	--static-quantized-quals 10 \
	--static-quantized-quals 20 \
	--static-quantized-quals 30 \
	-O $sample_id'.sorted.merged.MarkDups.BQSR.'$chunk_i'.bam'


}

chunk_process_s3(){
	i=$1
	chunk_i=$(basename $i ".bed")

	##GVCF
	$gatk_path/gatk --java-options "-Xmx4G" HaplotypeCaller \
	-R $assembly_ref \
	-ERC GVCF \
	-L $i \
	-OVI \
	-I $sample_id'.sorted.merged.MarkDups.BQSR.bam' \
	-O $sample_id'.sorted.merged.MarkDups.BQSR.'$chunk_i'.g.vcf.gz'

	tabix $sample_id'.sorted.merged.MarkDups.BQSR.'$chunk_i'.g.vcf.gz'
}

###split the genome into 10 chunks, and process the s1 for each chunk.
for j in $chunks_path/*.bed;do
	(
	chunk_process_s1 $j
	)&
done
wait

##combined the BQSR reports from each chunk
ls $sample_id'.sorted.merged.MarkDups.'*'.table' > $sample_id'.BQSR.reports.list'
#
$gatk_path/gatk --java-options "-Xmx6G" GatherBQSRReports \
-I $sample_id'.BQSR.reports.list' \
-O $sample_id'.BQSR.reports.table' && xargs -a $sample_id'.BQSR.reports.list' rm
#
rm $sample_id'.sorted.merged.bam'
rm $sample_id'.sorted.merged.bam.bai'
rm $sample_id'.sorted.merged.chunk_'*'.matrix'
rm $sample_id'.sorted.merged.chunk_'*'.bam'
rm $sample_id'.sorted.merged.chunk_'*'.bam.bai'
#
#
###apply BQSR with  combined table, and continue with s2 process
for j in $chunks_path/*.bed;do
	(
	chunk_process_s2 $j
	)&
done
wait

## merge the bam chunks
ls $sample_id'.sorted.merged.MarkDups.BQSR.chunk_'*'.bam' > chunks.bamlist
echo $sample_id'.unmapped.bam' >> chunks.bamlist

samtools merge -@ 10 -fb chunks.bamlist  $sample_id'.sorted.merged.MarkDups.BQSR.bam' &&
samtools index $sample_id'.sorted.merged.MarkDups.BQSR.bam' &&

rm chunks.bamlist
rm $sample_id'.sorted.merged.MarkDups.chunk_'*'.bam'
rm $sample_id'.sorted.merged.MarkDups.chunk_'*'.bam.bai'
rm $sample_id'.sorted.merged.MarkDups.BQSR.chunk_'*'.bam'
rm $sample_id'.sorted.merged.MarkDups.BQSR.chunk_'*'.bai'
rm $sample_id'.unmapped.bam'



##generate the g.vcf file for each chunk
for j in $chunks_path/*.bed;do
	(
	chunk_process_s3 $j
	)&
done
wait



##splite the gvcf for each chromosome, in order to merge them useing GatherVCFs

chunk_process_s4(){
	chunk_i='chunk_'$1
	tabix --list-chroms $sample_id'.sorted.merged.MarkDups.BQSR.'$chunk_i'.g.vcf.gz' > $chunk_i'.list'
	zcat $sample_id'.sorted.merged.MarkDups.BQSR.'$chunk_i'.g.vcf.gz' |grep "^#" > $chunk_i'.header'
	while IFS= read -r line;do
		tabix $sample_id'.sorted.merged.MarkDups.BQSR.'$chunk_i'.g.vcf.gz' $line |cat $chunk_i'.header' - |bgzip > $sample_id'.sorted.merged.MarkDups.BQSR.sc.'$line'.g.vcf.gz'
		tabix -p vcf $sample_id'.sorted.merged.MarkDups.BQSR.sc.'$line'.g.vcf.gz'
	done < $chunk_i'.list'
}

##no need to sort the contigs in chunk_11, just split the first 10 chunks

for j in {1..10};do
	(
	chunk_process_s4 $j
	)&
done
wait


##merge the g.vcf files

$gatk_path/gatk --java-options "-Xmx6G" GatherVcfs \
-I $sample_id'.sorted.merged.MarkDups.BQSR.sc.chr1.g.vcf.gz' \
-I $sample_id'.sorted.merged.MarkDups.BQSR.sc.chr2.g.vcf.gz' \
-I $sample_id'.sorted.merged.MarkDups.BQSR.sc.chr3.g.vcf.gz' \
-I $sample_id'.sorted.merged.MarkDups.BQSR.sc.chr4.g.vcf.gz' \
-I $sample_id'.sorted.merged.MarkDups.BQSR.sc.chr5.g.vcf.gz' \
-I $sample_id'.sorted.merged.MarkDups.BQSR.sc.chr6.g.vcf.gz' \
-I $sample_id'.sorted.merged.MarkDups.BQSR.sc.chr7.g.vcf.gz' \
-I $sample_id'.sorted.merged.MarkDups.BQSR.sc.chr8.g.vcf.gz' \
-I $sample_id'.sorted.merged.MarkDups.BQSR.sc.chr9.g.vcf.gz' \
-I $sample_id'.sorted.merged.MarkDups.BQSR.sc.chr10.g.vcf.gz' \
-I $sample_id'.sorted.merged.MarkDups.BQSR.sc.chr11.g.vcf.gz' \
-I $sample_id'.sorted.merged.MarkDups.BQSR.sc.chr12.g.vcf.gz' \
-I $sample_id'.sorted.merged.MarkDups.BQSR.sc.chr13.g.vcf.gz' \
-I $sample_id'.sorted.merged.MarkDups.BQSR.sc.chr14.g.vcf.gz' \
-I $sample_id'.sorted.merged.MarkDups.BQSR.sc.chr15.g.vcf.gz' \
-I $sample_id'.sorted.merged.MarkDups.BQSR.sc.chr16.g.vcf.gz' \
-I $sample_id'.sorted.merged.MarkDups.BQSR.sc.chr17.g.vcf.gz' \
-I $sample_id'.sorted.merged.MarkDups.BQSR.sc.chr18.g.vcf.gz' \
-I $sample_id'.sorted.merged.MarkDups.BQSR.sc.chr19.g.vcf.gz' \
-I $sample_id'.sorted.merged.MarkDups.BQSR.sc.chr20.g.vcf.gz' \
-I $sample_id'.sorted.merged.MarkDups.BQSR.sc.chr21.g.vcf.gz' \
-I $sample_id'.sorted.merged.MarkDups.BQSR.sc.chr22.g.vcf.gz' \
-I $sample_id'.sorted.merged.MarkDups.BQSR.sc.chr23.g.vcf.gz' \
-I $sample_id'.sorted.merged.MarkDups.BQSR.sc.chr24.g.vcf.gz' \
-I $sample_id'.sorted.merged.MarkDups.BQSR.sc.chr25.g.vcf.gz' \
-I $sample_id'.sorted.merged.MarkDups.BQSR.sc.chr26.g.vcf.gz' \
-I $sample_id'.sorted.merged.MarkDups.BQSR.sc.chr27.g.vcf.gz' \
-I $sample_id'.sorted.merged.MarkDups.BQSR.sc.chr28.g.vcf.gz' \
-I $sample_id'.sorted.merged.MarkDups.BQSR.sc.chr29.g.vcf.gz' \
-I $sample_id'.sorted.merged.MarkDups.BQSR.sc.chr30.g.vcf.gz' \
-I $sample_id'.sorted.merged.MarkDups.BQSR.sc.chr31.g.vcf.gz' \
-I $sample_id'.sorted.merged.MarkDups.BQSR.sc.chr32.g.vcf.gz' \
-I $sample_id'.sorted.merged.MarkDups.BQSR.sc.chr33.g.vcf.gz' \
-I $sample_id'.sorted.merged.MarkDups.BQSR.sc.chr34.g.vcf.gz' \
-I $sample_id'.sorted.merged.MarkDups.BQSR.sc.chr35.g.vcf.gz' \
-I $sample_id'.sorted.merged.MarkDups.BQSR.sc.chr36.g.vcf.gz' \
-I $sample_id'.sorted.merged.MarkDups.BQSR.sc.chr37.g.vcf.gz' \
-I $sample_id'.sorted.merged.MarkDups.BQSR.sc.chr38.g.vcf.gz' \
-I $sample_id'.sorted.merged.MarkDups.BQSR.sc.chrX.g.vcf.gz' \
-I $sample_id'.sorted.merged.MarkDups.BQSR.sc.chrY_NC_051844.1.g.vcf.gz' \
-I $sample_id'.sorted.merged.MarkDups.BQSR.sc.chrM.g.vcf.gz' \
-I $sample_id'.sorted.merged.MarkDups.BQSR.sc.chrY_unplaced_NW_024010443.1.g.vcf.gz' \
-I $sample_id'.sorted.merged.MarkDups.BQSR.sc.chrY_unplaced_NW_024010444.1.g.vcf.gz' \
-I $sample_id'.sorted.merged.MarkDups.BQSR.chunk_11.g.vcf.gz' \
-O $sample_id'.sorted.merged.MarkDups.BQSR.g.vcf.gz' &&

tabix -p vcf $sample_id'.sorted.merged.MarkDups.BQSR.g.vcf.gz'

##clean the files
rm chunk_*.list
rm chunk_*.header
rm $sample_id'.BQSR.reports.list'
rm $sample_id'.sorted.merged.MarkDups.BQSR.chunk_'*'.g.vcf.gz'
rm $sample_id'.sorted.merged.MarkDups.BQSR.chunk_'*'.g.vcf.gz.tbi'
rm $sample_id'.sorted.merged.MarkDups.BQSR.sc.chr'*'.g.vcf.gz'
rm $sample_id'.sorted.merged.MarkDups.BQSR.sc.chr'*'.g.vcf.gz.tbi'

python $run_stats --cram $sample_id'.sorted.merged.MarkDups.BQSR.bam' --ref $assembly_ref --gvcf $sample_id'.sorted.merged.MarkDups.BQSR.g.vcf.gz' --site $depth_sites

echo "Job is done"
date
