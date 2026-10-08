#!/bin/bash -l
#SBATCH -A $proj_id
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 15:00:00
#SBATCH -J finemapping_ldsc

## Compute annotation-specific LD scores
## https://github.com/bulik/ldsc/wiki/LD-Score-Estimation-Tutorial


dir_LDSC=dog-finemapping_polyfun/

source /sw/apps/conda/latest/rackham_stage/etc/profile.d/conda.sh
conda activate /home/kteng/.conda/envs/polyfun 

cd ${dir_LDSC}

#item146:
python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.1 \
    --sumstats CCD_250822/geno_files/ITEM146_sumstats_munged.parquet \
    --non-funct \
    --n 2419 \
    --chr 1 \
    --start 75042353 \
    --end 75250007 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM146_NO_PRIORS.1.200kb.gz

python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.15 \
    --sumstats CCD_250822/geno_files/ITEM146_sumstats_munged.parquet \
    --non-funct \
    --n 2419 \
    --chr 15 \
    --start 56998891 \
    --end 57198891 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM146_NO_PRIORS.15.200kb.gz


python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.3 \
    --sumstats CCD_250822/geno_files/ITEM146_sumstats_munged.parquet \
    --non-funct \
    --n 2419 \
    --chr 3 \
    --start 61192749 \
    --end 61394503 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM146_NO_PRIORS.3.200kb.gz

#item154

python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.4 \
    --sumstats CCD_250822/geno_files/ITEM154_sumstats_munged.parquet \
    --non-funct \
    --n 2418 \
    --chr 4 \
    --start 57387700 \
    --end 57752076 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM154_NO_PRIORS.4.360kb.gz

python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.4 \
    --sumstats CCD_250822/geno_files/ITEM154_sumstats_munged.parquet \
    --non-funct \
    --n 2418 \
    --chr 4 \
    --start 57487700 \
    --end 57652076 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM154_NO_PRIORS.4.160kb.gz


python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.15 \
    --sumstats CCD_250822/geno_files/ITEM154_sumstats_munged.parquet \
    --non-funct \
    --n 2418 \
    --chr 15 \
    --start 26426847 \
    --end 26643821 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM154_NO_PRIORS.15.220kb.gz

#########chr17 ########
python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.17 \
    --sumstats CCD_250822/geno_files/ITEM154_sumstats_munged.parquet \
    --non-funct \
    --n 2418 \
    --chr 17 \
    --start 4812399 \
    --end 5321576 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM154_NO_PRIORS.17.500kb.gz

python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.17 \
    --sumstats CCD_250822/geno_files/ITEM154_sumstats_munged.parquet \
    --non-funct \
    --n 2418 \
    --chr 17 \
    --start 4912399 \
    --end 5221576 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM154_NO_PRIORS.17.300kb.gz

#item153
python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.17 \
    --sumstats CCD_250822/geno_files/ITEM153_sumstats_munged.parquet \
    --non-funct \
    --n 2417 \
    --chr 17 \
    --start 4967580 \
    --end 5320702 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM153_NO_PRIORS.17.350kb.gz

python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.17 \
    --sumstats CCD_250822/geno_files/ITEM153_sumstats_munged.parquet \
    --non-funct \
    --n 2417 \
    --chr 17 \
    --start 5067580 \
    --end 5220702 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM153_NO_PRIORS.17.150kb.gz


python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.17 \
    --sumstats CCD_250822/geno_files/CCDF1_sumstats_munged.parquet \
    --non-funct \
    --n 2429 \
    --chr 17 \
    --start 5067510 \
    --end 5220702 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/CCDF1_NO_PRIORS.17.150kb.gz

python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.17 \
    --sumstats CCD_250822/geno_files/CCDF1_sumstats_munged.parquet \
    --non-funct \
    --n 2429 \
    --chr 17 \
    --start 4967510 \
    --end 5320702 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/CCDF1_NO_PRIORS.17.350kb.gz

########
#item154 cont.
python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.20 \
    --sumstats CCD_250822/geno_files/ITEM154_sumstats_munged.parquet \
    --non-funct \
    --n 2418 \
    --chr 20 \
    --start 11880789 \
    --end 12106620 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM154_NO_PRIORS.20.220kb.gz

python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.22 \
    --sumstats CCD_250822/geno_files/ITEM154_sumstats_munged.parquet \
    --non-funct \
    --n 2418 \
    --chr 22 \
    --start 25205335 \
    --end 25480179 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM154_NO_PRIORS.22_region1.274kb.gz

python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.22 \
    --sumstats CCD_250822/geno_files/ITEM154_sumstats_munged.parquet \
    --non-funct \
    --n 2418 \
    --chr 22 \
    --start 25365201 \
    --end 25855325 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM154_NO_PRIORS.22_region2.500kb.gz

python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.24 \
    --sumstats CCD_250822/geno_files/ITEM154_sumstats_munged.parquet \
    --non-funct \
    --n 2418 \
    --chr 24 \
    --start 28855671 \
    --end 29060651 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM154_NO_PRIORS.24.200kb.gz


python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.24 \
    --sumstats CCD_250822/geno_files/ITEM154_sumstats_munged.parquet \
    --non-funct \
    --n 2418 \
    --chr 24 \
    --start 28955671 \
    --end 28960651 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM154_NO_PRIORS.24.5kb.gz