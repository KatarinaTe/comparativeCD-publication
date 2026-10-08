#!/bin/bash -l
#SBATCH -A snic2022-5-416
#SBATCH -p core -n 1
#SBATCH -J GlimpseChunk
#SBATCH -t 1:00:00
#SBATCH -o /path/to/chunk.output
#SBATCH -e /path/to/Z_ErrorAndOut/chunk.error
#SBATCH --mail-user 
#SBATCH --mail-type=ALL

module load bioinfo-tools

ImputationRefFolder='/path/to/1_phased-imputation-panel'
WorkFolder='/path/to/GLIMPSEv1.1/GLIMPSE_chr_chunks'
ScriptsFolder='/path/to/0_Scripts'

cd /proj/snic2022-6-229/nobackup/private/Jennifer/GLIMPSEv1.1

./GLIMPSE_chunk_static --input $ImputationRefFolder/AutoAndXPAR.Dog10K.Phased.bcf --region chr1 --output $WorkFolder/chunks_chr1.txt
./GLIMPSE_chunk_static --input $ImputationRefFolder/AutoAndXPAR.Dog10K.Phased.bcf --region chr2 --output $WorkFolder/chunks_chr2.txt
./GLIMPSE_chunk_static --input $ImputationRefFolder/AutoAndXPAR.Dog10K.Phased.bcf --region chr3 --output $WorkFolder/chunks_chr3.txt
./GLIMPSE_chunk_static --input $ImputationRefFolder/AutoAndXPAR.Dog10K.Phased.bcf --region chr4 --output $WorkFolder/chunks_chr4.txt
./GLIMPSE_chunk_static --input $ImputationRefFolder/AutoAndXPAR.Dog10K.Phased.bcf --region chr5 --output $WorkFolder/chunks_chr5.txt
./GLIMPSE_chunk_static --input $ImputationRefFolder/AutoAndXPAR.Dog10K.Phased.bcf --region chr6 --output $WorkFolder/chunks_chr6.txt
./GLIMPSE_chunk_static --input $ImputationRefFolder/AutoAndXPAR.Dog10K.Phased.bcf --region chr7 --output $WorkFolder/chunks_chr7.txt
./GLIMPSE_chunk_static --input $ImputationRefFolder/AutoAndXPAR.Dog10K.Phased.bcf --region chr8 --output $WorkFolder/chunks_chr8.txt
./GLIMPSE_chunk_static --input $ImputationRefFolder/AutoAndXPAR.Dog10K.Phased.bcf --region chr9 --output $WorkFolder/chunks_chr9.txt
./GLIMPSE_chunk_static --input $ImputationRefFolder/AutoAndXPAR.Dog10K.Phased.bcf --region chr10 --output $WorkFolder/chunks_chr10.txt
./GLIMPSE_chunk_static --input $ImputationRefFolder/AutoAndXPAR.Dog10K.Phased.bcf --region chr11 --output $WorkFolder/chunks_chr11.txt
./GLIMPSE_chunk_static --input $ImputationRefFolder/AutoAndXPAR.Dog10K.Phased.bcf --region chr12 --output $WorkFolder/chunks_chr12.txt
./GLIMPSE_chunk_static --input $ImputationRefFolder/AutoAndXPAR.Dog10K.Phased.bcf --region chr13 --output $WorkFolder/chunks_chr13.txt
./GLIMPSE_chunk_static --input $ImputationRefFolder/AutoAndXPAR.Dog10K.Phased.bcf --region chr14 --output $WorkFolder/chunks_chr14.txt
./GLIMPSE_chunk_static --input $ImputationRefFolder/AutoAndXPAR.Dog10K.Phased.bcf --region chr15 --output $WorkFolder/chunks_chr15.txt
./GLIMPSE_chunk_static --input $ImputationRefFolder/AutoAndXPAR.Dog10K.Phased.bcf --region chr16 --output $WorkFolder/chunks_chr16.txt
./GLIMPSE_chunk_static --input $ImputationRefFolder/AutoAndXPAR.Dog10K.Phased.bcf --region chr17 --output $WorkFolder/chunks_chr17.txt
./GLIMPSE_chunk_static --input $ImputationRefFolder/AutoAndXPAR.Dog10K.Phased.bcf --region chr18 --output $WorkFolder/chunks_chr18.txt
./GLIMPSE_chunk_static --input $ImputationRefFolder/AutoAndXPAR.Dog10K.Phased.bcf --region chr19 --output $WorkFolder/chunks_chr19.txt
./GLIMPSE_chunk_static --input $ImputationRefFolder/AutoAndXPAR.Dog10K.Phased.bcf --region chr20 --output $WorkFolder/chunks_chr20.txt
./GLIMPSE_chunk_static --input $ImputationRefFolder/AutoAndXPAR.Dog10K.Phased.bcf --region chr21 --output $WorkFolder/chunks_chr21.txt
./GLIMPSE_chunk_static --input $ImputationRefFolder/AutoAndXPAR.Dog10K.Phased.bcf --region chr22 --output $WorkFolder/chunks_chr22.txt
./GLIMPSE_chunk_static --input $ImputationRefFolder/AutoAndXPAR.Dog10K.Phased.bcf --region chr23 --output $WorkFolder/chunks_chr23.txt
./GLIMPSE_chunk_static --input $ImputationRefFolder/AutoAndXPAR.Dog10K.Phased.bcf --region chr24 --output $WorkFolder/chunks_chr24.txt
./GLIMPSE_chunk_static --input $ImputationRefFolder/AutoAndXPAR.Dog10K.Phased.bcf --region chr25 --output $WorkFolder/chunks_chr25.txt
./GLIMPSE_chunk_static --input $ImputationRefFolder/AutoAndXPAR.Dog10K.Phased.bcf --region chr26 --output $WorkFolder/chunks_chr26.txt
./GLIMPSE_chunk_static --input $ImputationRefFolder/AutoAndXPAR.Dog10K.Phased.bcf --region chr27 --output $WorkFolder/chunks_chr27.txt
./GLIMPSE_chunk_static --input $ImputationRefFolder/AutoAndXPAR.Dog10K.Phased.bcf --region chr28 --output $WorkFolder/chunks_chr28.txt
./GLIMPSE_chunk_static --input $ImputationRefFolder/AutoAndXPAR.Dog10K.Phased.bcf --region chr29 --output $WorkFolder/chunks_chr29.txt
./GLIMPSE_chunk_static --input $ImputationRefFolder/AutoAndXPAR.Dog10K.Phased.bcf --region chr30 --output $WorkFolder/chunks_chr30.txt
./GLIMPSE_chunk_static --input $ImputationRefFolder/AutoAndXPAR.Dog10K.Phased.bcf --region chr31 --output $WorkFolder/chunks_chr31.txt
./GLIMPSE_chunk_static --input $ImputationRefFolder/AutoAndXPAR.Dog10K.Phased.bcf --region chr32 --output $WorkFolder/chunks_chr32.txt
./GLIMPSE_chunk_static --input $ImputationRefFolder/AutoAndXPAR.Dog10K.Phased.bcf --region chr33 --output $WorkFolder/chunks_chr33.txt
./GLIMPSE_chunk_static --input $ImputationRefFolder/AutoAndXPAR.Dog10K.Phased.bcf --region chr34 --output $WorkFolder/chunks_chr34.txt
./GLIMPSE_chunk_static --input $ImputationRefFolder/AutoAndXPAR.Dog10K.Phased.bcf --region chr35 --output $WorkFolder/chunks_chr35.txt
./GLIMPSE_chunk_static --input $ImputationRefFolder/AutoAndXPAR.Dog10K.Phased.bcf --region chr36 --output $WorkFolder/chunks_chr36.txt
./GLIMPSE_chunk_static --input $ImputationRefFolder/AutoAndXPAR.Dog10K.Phased.bcf --region chr37 --output $WorkFolder/chunks_chr37.txt
./GLIMPSE_chunk_static --input $ImputationRefFolder/AutoAndXPAR.Dog10K.Phased.bcf --region chr38 --output $WorkFolder/chunks_chr38.txt
