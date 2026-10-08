#!/bin/bash -l
#SBATCH -A $proj_id
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 15:00:00
#SBATCH -J finemapping_ldsc

## Compute annotation-specific LD scores
## https://github.com/bulik/ldsc/wiki/LD-Score-Estimation-Tutorial


dir_LDSC=/dog-finemapping_polyfun/

source /sw/apps/conda/latest/rackham_stage/etc/profile.d/conda.sh
conda activate /home/kteng/.conda/envs/polyfun 

cd ${dir_LDSC}

#item149:

python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.1 \
    --sumstats CCD_250822/geno_files/ITEM149_sumstats_munged.parquet \
    --non-funct \
    --n 2422 \
    --chr 1 \
    --start 117371833 \
    --end 117853446 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM149_NO_PRIORS.1.480kb.gz

python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.11 \
    --sumstats CCD_250822/geno_files/ITEM149_sumstats_munged.parquet \
    --non-funct \
    --n 2422 \
    --chr 11 \
    --start 50174733 \
    --end 50374733 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM149_NO_PRIORS.11.200kb.gz

python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.13 \
    --sumstats CCD_250822/geno_files/ITEM149_sumstats_munged.parquet \
    --non-funct \
    --n 2422 \
    --chr 13 \
    --start 41931759 \
    --end 42169584 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM149_NO_PRIORS.13.240kb.gz

##narrowed:
python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.13 \
    --sumstats CCD_250822/geno_files/ITEM149_sumstats_munged.parquet \
    --non-funct \
    --n 2422 \
    --chr 13 \
    --start 42031759 \
    --end 42069584 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM149_NO_PRIORS.13.40kb.gz



python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.26 \
    --sumstats CCD_250822/geno_files/ITEM149_sumstats_munged.parquet \
    --non-funct \
    --n 2422 \
    --chr 26 \
    --start 29115127 \
    --end 29316918 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM149_NO_PRIORS.26.200kb.gz 


python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.29 \
    --sumstats CCD_250822/geno_files/ITEM149_sumstats_munged.parquet \
    --non-funct \
    --n 2422 \
    --chr 29 \
    --start 36799088 \
    --end 36999341 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM149_NO_PRIORS.29.200kb.gz

python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.29 \
    --sumstats CCD_250822/geno_files/ITEM149_sumstats_munged.parquet \
    --non-funct \
    --n 2422 \
    --chr 29 \
    --start 36899088 \
    --end 36899341 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM149_NO_PRIORS.29.0.3kb.gz




python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.38 \
    --sumstats CCD_250822/geno_files/ITEM149_sumstats_munged.parquet \
    --non-funct \
    --n 2422 \
    --chr 38 \
    --start 19152704 \
    --end 19352843 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM149_NO_PRIORS.38.200kb.gz

python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.6 \
    --sumstats CCD_250822/geno_files/ITEM149_sumstats_munged.parquet \
    --non-funct \
    --n 2422 \
    --chr 6 \
    --start 23030194 \
    --end 23238102 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM149_NO_PRIORS.6.200kb.gz

#item145:
python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.2 \
    --sumstats CCD_250822/geno_files/ITEM145_sumstats_munged.parquet \
    --non-funct \
    --n 2418 \
    --chr 2 \
    --start 57776761 \
    --end 58034887 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM145_NO_PRIORS.2_region1.260kb.gz 

python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.2 \
    --sumstats CCD_250822/geno_files/ITEM145_sumstats_munged.parquet \
    --non-funct \
    --n 2418 \
    --chr 2 \
    --start 57876761 \
    --end 57934887 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM145_NO_PRIORS.2_region1.60kb.gz 

python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.2 \
    --sumstats CCD_250822/geno_files/ITEM145_sumstats_munged.parquet \
    --non-funct \
    --n 2418 \
    --chr 2 \
    --start 57790238 \
    --end 58181791 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM145_NO_PRIORS.2_region2.400kb.gz 

python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.2 \
    --sumstats CCD_250822/geno_files/ITEM145_sumstats_munged.parquet \
    --non-funct \
    --n 2418 \
    --chr 2 \
    --start 57890238 \
    --end 58081791 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM145_NO_PRIORS.2_region2.200kb.gz 


python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.4 \
    --sumstats CCD_250822/geno_files/ITEM145_sumstats_munged.parquet \
    --non-funct \
    --n 2418 \
    --chr 4 \
    --start 76977659 \
    --end 77177659 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM145_NO_PRIORS.4.200kb.gz 

python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.7 \
    --sumstats CCD_250822/geno_files/ITEM145_sumstats_munged.parquet \
    --non-funct \
    --n 2418 \
    --chr 7 \
    --start 6499147 \
    --end 6708947 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM145_NO_PRIORS.7_region1.210kb.gz

python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.7 \
    --sumstats CCD_250822/geno_files/ITEM145_sumstats_munged.parquet \
    --non-funct \
    --n 2418 \
    --chr 7 \
    --start 6599147 \
    --end 6608947 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM145_NO_PRIORS.7_region1.10kb.gz

python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.7 \
    --sumstats CCD_250822/geno_files/ITEM145_sumstats_munged.parquet \
    --non-funct \
    --n 2418 \
    --chr 7 \
    --start 6508991 \
    --end 6708991 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM145_NO_PRIORS.7_region2.200kb.gz

python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.7 \
    --sumstats CCD_250822/geno_files/ITEM145_sumstats_munged.parquet \
    --non-funct \
    --n 2418 \
    --chr 7 \
    --start 60997669 \
    --end 61197674 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM145_NO_PRIORS.7_region3.200kb.gz

python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.11 \
    --sumstats CCD_250822/geno_files/ITEM145_sumstats_munged.parquet \
    --non-funct \
    --n 2418 \
    --chr 11 \
    --start 72828242 \
    --end 73028242 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM145_NO_PRIORS.11.200kb.gz 

python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.12 \
    --sumstats CCD_250822/geno_files/ITEM145_sumstats_munged.parquet \
    --non-funct \
    --n 2418 \
    --chr 12 \
    --start 11278620 \
    --end 11587660 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM145_NO_PRIORS.12_region1.310kb.gz 


python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.12 \
    --sumstats CCD_250822/geno_files/ITEM145_sumstats_munged.parquet \
    --non-funct \
    --n 2418 \
    --chr 12 \
    --start 33619438 \
    --end 33868329 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM145_NO_PRIORS.12_region2.250kb.gz 


python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.13 \
    --sumstats CCD_250822/geno_files/ITEM145_sumstats_munged.parquet \
    --non-funct \
    --n 2418 \
    --chr 13 \
    --start 9579944 \
    --end 9918047 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM145_NO_PRIORS.13.300kb.gz

python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.32 \
    --sumstats CCD_250822/geno_files/ITEM145_sumstats_munged.parquet \
    --non-funct \
    --n 2418 \
    --chr 32 \
    --start 13519324 \
    --end 13743397 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM145_NO_PRIORS.32.13.220kb.gz


python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.35 \
    --sumstats CCD_250822/geno_files/ITEM145_sumstats_munged.parquet \
    --non-funct \
    --n 2418 \
    --chr 35 \
    --start 8295148 \
    --end 8499503 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM145_NO_PRIORS.35_region1.200kb.gz

python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.35 \
    --sumstats CCD_250822/geno_files/ITEM145_sumstats_munged.parquet \
    --non-funct \
    --n 2418 \
    --chr 35 \
    --start 8300366 \
    --end 8502086 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM145_NO_PRIORS.35_region2.200kb.gz

