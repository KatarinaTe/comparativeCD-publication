#SBATCH -A $proj_id
#SBATCH -p core
#SBATCH -J ImputeAndPhaseChrsByChunks
#SBATCH -t 24:00:00

module load bioinfo-tools
module load bcftools
module load samtools

input_file=$1


ScriptsFolder='/path/to/0_Scripts/ChunkList.Chr3-10' #OBS, change script for each ChunkList:
#ScriptsFolder='/path/to/0_Scripts/ChunkList.Chr1'
#ScriptsFolder='/path/to/0_Scripts/ChunkList.Chr2'
#ScriptsFolder='/path/to/0_Scripts/ChunkList.Chr11-20'
#ScriptsFolder='/path/to/0_Scripts/ChunkList.Chr21-30'
#ScriptsFolder='/path/to/0_Scripts/ChunkList.Chr31-37'
#ScriptsFolder='/path/to/0_Scripts/ChunkList.Chr38' #I use one folder for Chr38 even if that is the smallest, to try that everything works.


SAMPLEVCFFolder='/path/to/3_GenotypeLikelihood' # GenotypeLikihood output
Reference='/path/to/1_phased-imputation-panel/AutoAndXPAR.Dog10K.Phased.bcf' #original phased panel, input for step 3 when extracting SNPs 

PHASEOUT='/path/to/4_ImputeAndPhase/ChunkList.Chr3-10' #OBS, change for each ChunkList:
#PHASEOUT='/path/to/4_ImputeAndPhase/ChunkList.Chr1'
#PHASEOUT='/path/to/4_ImputeAndPhase/ChunkList.Chr2'
#PHASEOUT='/path/to/4_ImputeAndPhase/ChunkList.Chr11-20'
#PHASEOUT='/path/to/4_ImputeAndPhase/ChunkList.Chr21-30'
#PHASEOUT='/path/to/4_ImputeAndPhase/ChunkList.Chr31-37'
#PHASEOUT='/path/to/4_ImputeAndPhase/ChunkList.Chr38' 

CHUNKS='/path/to/GLIMPSEv1.1/GLIMPSE_chr_chunks'

# To get lines that are needed for input -> in /path/to/0_Scripts/ChunkList.Chr3-10/
# split -l 1 /path/to/GLIMPSEv1.1/GLIMPSE_chr_chunks/chunks_chr3.txt chr3
# split -l 1 /path/to/GLIMPSEv1.1/GLIMPSE_chr_chunks/chunks_chr4.txt chr4
#...until chr10
#then do the same for the other ChunkList folders.


# within Chunks folder where just made lines, i.e., /path/to/0_Scripts/ChunkList.Chr3-10/
# Script Usage -> for i in chr*;do echo $i; sbatch ImputeAndPhaseChrsByChuncks_Chr3-10.sh $i; done
#then do the same for the other ChunkList folders.


cd /path/to/GLIMPSEv1.1/

#OBS!! remember to change naming of samples, here DA.LowPass.vcf.gz

while IFS="" read -r LINE || [ -n "$LINE" ]; 
do   
        printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
        CHR=$(echo $LINE | cut -d" " -f2)
        IRG=$(echo $LINE | cut -d" " -f3)
        ORG=$(echo $LINE | cut -d" " -f4)
        OUT=${PHASEOUT}/DA.${CHR}.imputed.${ID}.bcf
        ./GLIMPSE_phase_static --input ${SAMPLEVCFFolder}/${CHR}.GL.DA_LowPass.vcf.gz --reference ${Reference} --input-region ${IRG} --output-region ${ORG} --output ${OUT}
        bcftools index -f ${OUT}
done < $ScriptsFolder/$input_file