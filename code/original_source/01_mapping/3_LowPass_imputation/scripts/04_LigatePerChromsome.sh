#!/bin/bash -l
#SBATCH -A $proj_id
#SBATCH -p core -n 1
#SBATCH -J LigatePerChromsome
#SBATCH -t 10:00:00
#SBATCH -o /path/to/Z_ErrorAndOut/LigatePerChromsome.output
#SBATCH -e /path/to/Z_ErrorAndOut/LigatePerChromsome.error
#SBATCH --mail-user $user_email
#SBATCH --mail-type=ALL

module load bioinfo-tools
module load bcftools
module load samtools

Reference='/path/to/UU_Cfam_GSD_1.0_ROSY.fa'
PhaseChunk='/path/to/4_ImputeAndPhase/ChunkList.Chr1'
#PhaseChunk='/path/to/4_ImputeAndPhase/ChunkList.Chr2'
#PhaseChunk='/path/to/4_ImputeAndPhase/ChunkList.Chr3-10'
#PhaseChunk='/path/to/4_ImputeAndPhase/ChunkList.Chr11-20'
#PhaseChunk='/path/to/4_ImputeAndPhase/ChunkList.Chr21-30'
#PhaseChunk='/path/to/4_ImputeAndPhase/ChunkList.Chr31-37'
#PhaseChunk='/path/to/4_ImputeAndPhase/ChunkList.Chr38' #only chr38 to use as trial run

LigateCHRs='/path/to/5_LigatedImputedCHRs'
ScriptsFolder='/path/to/0_Scripts'
Glimpse='/path/to/GLIMPSEv1.1'

cd $PhaseChunk

#OBS!! change chr$ according to the ChunkList!!!

$Glimpse/GLIMPSE_ligate_static --input list.chr1.txt --output $LigateCHRs/DA.chr1.merged.bcf
#$Glimpse/GLIMPSE_ligate_static --input list.chr2.txt --output $LigateCHRs/DA.chr2.merged.bcf
#$Glimpse/GLIMPSE_ligate_static --input list.chr3.txt --output $LigateCHRs/DA.chr3.merged.bcf
#$Glimpse/GLIMPSE_ligate_static --input list.chr4.txt --output $LigateCHRs/DA.chr4.merged.bcf
#etc

# The script above makes an index for the bcf
 
bcftools convert -O b $LigateCHRs/DA.chr1.merged.bcf -o $LigateCHRs/DA.chr1.merged.vcf
bgzip -c $LigateCHRs/DA.chr1.merged.vcf > $LigateCHRs/DA.chr1.merged.vcf.gz
bcftools index -f $LigateCHRs/DA.chr1.merged.vcf.gz

#bcftools convert -O b $LigateCHRs/DA.chr2.merged.bcf -o $LigateCHRs/DA.chr2.merged.vcf
#bgzip -c $LigateCHRs/DA.chr2.merged.vcf > $LigateCHRs/DA.chr2.merged.vcf.gz
#bcftools index -f $LigateCHRs/DA.chr2.merged.vcf.gz

#bcftools convert -O b $LigateCHRs/DA.chr3.merged.bcf -o $LigateCHRs/DA.chr3.merged.vcf
#bgzip -c $LigateCHRs/DA.chr3.merged.vcf > $LigateCHRs/DA.chr3.merged.vcf.gz
#bcftools index -f $LigateCHRs/DA.chr3.merged.vcf.gz

#etc.