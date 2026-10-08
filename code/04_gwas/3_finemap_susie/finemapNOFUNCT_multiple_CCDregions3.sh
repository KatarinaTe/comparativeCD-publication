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


###continue from here:
#item155
python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.14 \
    --sumstats CCD_250822/geno_files/ITEM155_sumstats_munged.parquet \
    --non-funct \
    --n 2322 \
    --chr 14 \
    --start 58094237 \
    --end 58296796 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM155_NO_PRIORS.14.200kb.gz

python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.15 \
    --sumstats CCD_250822/geno_files/ITEM155_sumstats_munged.parquet \
    --non-funct \
    --n 2322 \
    --chr 15 \
    --start 40459377 \
    --end 40684487 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM155_NO_PRIORS.15.200kb.gz

python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.15 \
    --sumstats CCD_250822/geno_files/ITEM155_sumstats_munged.parquet \
    --non-funct \
    --n 2322 \
    --chr 15 \
    --start 40559377 \
    --end 40584487 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM155_NO_PRIORS.15.25kb.gz

python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.18 \
    --sumstats CCD_250822/geno_files/ITEM155_sumstats_munged.parquet \
    --non-funct \
    --n 2322 \
    --chr 18 \
    --start 52394924 \
    --end 52666653 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM155_NO_PRIORS.18.270kb.gz

python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.19 \
    --sumstats CCD_250822/geno_files/ITEM155_sumstats_munged.parquet \
    --non-funct \
    --n 2322 \
    --chr 19 \
    --start 53196143 \
    --end 53396143 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM155_NO_PRIORS.19.200kb.gz

python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.20 \
    --sumstats CCD_250822/geno_files/ITEM155_sumstats_munged.parquet \
    --non-funct \
    --n 2322 \
    --chr 20 \
    --start 7518972 \
    --end 7720413 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM155_NO_PRIORS.20.200kb.gz


python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.22 \
    --sumstats CCD_250822/geno_files/ITEM155_sumstats_munged.parquet \
    --non-funct \
    --n 2322 \
    --chr 22 \
    --start 19052988 \
    --end 19269451 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM155_NO_PRIORS.22.220kb.gz

python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.27 \
    --sumstats CCD_250822/geno_files/ITEM155_sumstats_munged.parquet \
    --non-funct \
    --n 2322 \
    --chr 27 \
    --start 43795819 \
    --end 43995819 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM155_NO_PRIORS.27.200kb.gz

python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.34 \
    --sumstats CCD_250822/geno_files/ITEM155_sumstats_munged.parquet \
    --non-funct \
    --n 2322 \
    --chr 34 \
    --start 40500334 \
    --end 40700334 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM155_NO_PRIORS.34.200kb.gz


python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.37 \
    --sumstats CCD_250822/geno_files/ITEM155_sumstats_munged.parquet \
    --non-funct \
    --n 2322 \
    --chr 37 \
    --start 25940092 \
    --end 26140256 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM155_NO_PRIORS.37.200kb.gz

python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.37 \
    --sumstats CCD_250822/geno_files/ITEM155_sumstats_munged.parquet \
    --non-funct \
    --n 2322 \
    --chr 37 \
    --start 26040092 \
    --end 26040256 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM155_NO_PRIORS.37.0.2kb.gz


#item7
python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.2 \
    --sumstats CCD_250822/geno_files/ITEM7_sumstats_munged.parquet \
    --non-funct \
    --n 2555 \
    --chr 2 \
    --start 31264801 \
    --end 31464801 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM7_NO_PRIORS.2.200kb.gz

python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.3 \
    --sumstats CCD_250822/geno_files/ITEM7_sumstats_munged.parquet \
    --non-funct \
    --n 2555 \
    --chr 3 \
    --start 89805630 \
    --end 90169186 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM7_NO_PRIORS.3.360kb.gz

python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.5 \
    --sumstats CCD_250822/geno_files/ITEM7_sumstats_munged.parquet \
    --non-funct \
    --n 2555 \
    --chr 5 \
    --start 79027032 \
    --end 79269269 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM7_NO_PRIORS.5.240kb.gz

python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.19 \
    --sumstats CCD_250822/geno_files/ITEM7_sumstats_munged.parquet \
    --non-funct \
    --n 2555 \
    --chr 19 \
    --start 47343491 \
    --end 47551269 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM7_NO_PRIORS.19.200kb.gz

python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.37 \
    --sumstats CCD_250822/geno_files/ITEM7_sumstats_munged.parquet \
    --non-funct \
    --n 2555 \
    --chr 37 \
    --start 26040092 \
    --end 26040256 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM7_NO_PRIORS.37.0.2kb.gz

python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.27 \
    --sumstats CCD_250822/geno_files/ITEM7_sumstats_munged.parquet \
    --non-funct \
    --n 2555 \
    --chr 27 \
    --start 18871774 \
    --end 19161955 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM7_NO_PRIORS.27.290kb.gz


python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.32 \
    --sumstats CCD_250822/geno_files/ITEM7_sumstats_munged.parquet \
    --non-funct \
    --n 2555 \
    --chr 32 \
    --start 36150418 \
    --end 36350505 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM7_NO_PRIORS.32.200kb.gz

python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.32 \
    --sumstats CCD_250822/geno_files/ITEM7_sumstats_munged.parquet \
    --non-funct \
    --n 2555 \
    --chr 32 \
    --start 36250418 \
    --end 36250505 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM7_NO_PRIORS.32.0.1kb.gz

python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.34 \
    --sumstats CCD_250822/geno_files/ITEM7_sumstats_munged.parquet \
    --non-funct \
    --n 2555 \
    --chr 34 \
    --start 21197914 \
    --end 21397914 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM7_NO_PRIORS.34.200kb.gz 

python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.35 \
    --sumstats CCD_250822/geno_files/ITEM7_sumstats_munged.parquet \
    --non-funct \
    --n 2555 \
    --chr 35 \
    --start 3541072 \
    --end 3902576 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM7_NO_PRIORS.35.region1.360kb.gz

python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.35 \
    --sumstats CCD_250822/geno_files/ITEM7_sumstats_munged.parquet \
    --non-funct \
    --n 2555 \
    --chr 35 \
    --start 3641072 \
    --end 3802576 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM7_NO_PRIORS.35.region1.160kb.gz

python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.35 \
    --sumstats CCD_250822/geno_files/ITEM7_sumstats_munged.parquet \
    --non-funct \
    --n 2555 \
    --chr 35 \
    --start 3761821 \
    --end 4037722 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM7_NO_PRIORS.35.region2.280kb.gz

python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.35 \
    --sumstats CCD_250822/geno_files/ITEM7_sumstats_munged.parquet \
    --non-funct \
    --n 2555 \
    --chr 35 \
    --start 3861821 \
    --end 3937722 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM7_NO_PRIORS.35.region2.80kb.gz