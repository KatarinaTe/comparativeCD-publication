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
#item93

python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.7 \
    --sumstats CCD_250822/geno_files/ITEM93_sumstats_munged.parquet \
    --non-funct \
    --n 2494 \
    --chr 7 \
    --start 51854674 \
    --end 52094198 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM93_NO_PRIORS.7.240kb.gz

python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.17 \
    --sumstats CCD_250822/geno_files/ITEM93_sumstats_munged.parquet \
    --non-funct \
    --n 2494 \
    --chr 17 \
    --start 30392239 \
    --end 30607446 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM93_NO_PRIORS.17.region1.220kb.gz

python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.17 \
    --sumstats CCD_250822/geno_files/ITEM93_sumstats_munged.parquet \
    --non-funct \
    --n 2494 \
    --chr 17 \
    --start 30492239 \
    --end 30507446 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM93_NO_PRIORS.17.region1.20kb.gz

python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.17 \
    --sumstats CCD_250822/geno_files/ITEM93_sumstats_munged.parquet \
    --non-funct \
    --n 2494 \
    --chr 17 \
    --start 31897137 \
    --end 32097137 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM93_NO_PRIORS.17.region2.200kb.gz


python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.19 \
    --sumstats CCD_250822/geno_files/ITEM93_sumstats_munged.parquet \
    --non-funct \
    --n 2494 \
    --chr 19 \
    --start 9799439 \
    --end 10005259 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM93_NO_PRIORS.19.200kb.gz

python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.19 \
    --sumstats CCD_250822/geno_files/ITEM93_sumstats_munged.parquet \
    --non-funct \
    --n 2494 \
    --chr 19 \
    --start 9899439 \
    --end 9905259 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM93_NO_PRIORS.19.6kb.gz

python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.30 \
    --sumstats CCD_250822/geno_files/ITEM93_sumstats_munged.parquet \
    --non-funct \
    --n 2494 \
    --chr 30 \
    --start 11709732 \
    --end 11919385 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM93_NO_PRIORS.30.210kb.gz


python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.31 \
    --sumstats CCD_250822/geno_files/ITEM93_sumstats_munged.parquet \
    --non-funct \
    --n 2494 \
    --chr 31 \
    --start 15100674 \
    --end 15310134 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM93_NO_PRIORS.31.210kb.gz


#item95
python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.5 \
    --sumstats CCD_250822/geno_files/ITEM95_sumstats_munged.parquet \
    --non-funct \
    --n 2498 \
    --chr 5 \
    --start 63030588 \
    --end 63248313 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM95_NO_PRIORS.5.210kb.gz

python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.5 \
    --sumstats CCD_250822/geno_files/ITEM95_sumstats_munged.parquet \
    --non-funct \
    --n 2498 \
    --chr 5 \
    --start 63130588 \
    --end 63148313 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM95_NO_PRIORS.5.17kb.gz


#item153
python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.1 \
    --sumstats CCD_250822/geno_files/ITEM153_sumstats_munged.parquet \
    --non-funct \
    --n 2417 \
    --chr 1 \
    --start 120975788 \
    --end 121175788 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM153_NO_PRIORS.1.200kb.gz


python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.6 \
    --sumstats CCD_250822/geno_files/ITEM153_sumstats_munged.parquet \
    --non-funct \
    --n 2417 \
    --chr 6 \
    --start 19178651 \
    --end 19398424 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM153_NO_PRIORS.6.220kb.gz

python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.6 \
    --sumstats CCD_250822/geno_files/ITEM153_sumstats_munged.parquet \
    --non-funct \
    --n 2417 \
    --chr 6 \
    --start 19278651 \
    --end 19298424 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM153_NO_PRIORS.6.20kb.gz

python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.15 \
    --sumstats CCD_250822/geno_files/ITEM153_sumstats_munged.parquet \
    --non-funct \
    --n 2417 \
    --chr 15 \
    --start 26426847 \
    --end 26643821 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM153_NO_PRIORS.15.220kb.gz 


python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.17 \
    --sumstats CCD_250822/geno_files/ITEM153_sumstats_munged.parquet \
    --non-funct \
    --n 2417 \
    --chr 17 \
    --start 13043473 \
    --end 13244342 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM153_NO_PRIORS.17.201kb.gz

python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.17 \
    --sumstats CCD_250822/geno_files/ITEM153_sumstats_munged.parquet \
    --non-funct \
    --n 2417 \
    --chr 17 \
    --start 4655003 \
    --end 5084803 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM153_NO_PRIORS.17.430kb.gz



python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.20 \
    --sumstats CCD_250822/geno_files/ITEM153_sumstats_munged.parquet \
    --non-funct \
    --n 2417 \
    --chr 20 \
    --start 11880789 \
    --end 12127264 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM153_NO_PRIORS.20.250kb.gz

python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.22 \
    --sumstats CCD_250822/geno_files/ITEM153_sumstats_munged.parquet \
    --non-funct \
    --n 2417 \
    --chr 22 \
    --start 25209643 \
    --end 25480179 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM153_NO_PRIORS.22.region1.270kb.gz

python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.22 \
    --sumstats CCD_250822/geno_files/ITEM153_sumstats_munged.parquet \
    --non-funct \
    --n 2417 \
    --chr 22 \
    --start 56147712 \
    --end 56348310 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM153_NO_PRIORS.22.region2.201kb.gz

python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.26 \
    --sumstats CCD_250822/geno_files/ITEM153_sumstats_munged.parquet \
    --non-funct \
    --n 2417 \
    --chr 26 \
    --start 32866691 \
    --end 33069818 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM153_NO_PRIORS.26.203kb.gz


#item147
python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.8 \
    --sumstats CCD_250822/geno_files/ITEM147_sumstats_munged.parquet \
    --non-funct \
    --n 2414 \
    --chr 8 \
    --start 56421200 \
    --end 56625211 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM147_NO_PRIORS.8.204kb.gz

#item148
python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.21 \
    --sumstats CCD_250822/geno_files/ITEM148_sumstats_munged.parquet \
    --non-funct \
    --n 2414 \
    --chr 21 \
    --start 15231294 \
    --end 15451077 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM148_NO_PRIORS.21.220kb.gz

python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.27 \
    --sumstats CCD_250822/geno_files/ITEM148_sumstats_munged.parquet \
    --non-funct \
    --n 2414 \
    --chr 27 \
    --start 41559417 \
    --end 41760966 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM148_NO_PRIORS.27.201kb.gz

python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.27 \
    --sumstats CCD_250822/geno_files/ITEM148_sumstats_munged.parquet \
    --non-funct \
    --n 2414 \
    --chr 27 \
    --start 41659417 \
    --end 41660966 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM148_NO_PRIORS.27.1kb.gz



##CCDF1
python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.12 \
    --sumstats CCD_250822/geno_files/CCDF1_sumstats_munged.parquet \
    --non-funct \
    --n 2429 \
    --chr 12 \
    --start 3384409 \
    --end 3585864 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/CCDF1_NO_PRIORS.12.201kb.gz

python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.15 \
    --sumstats CCD_250822/geno_files/CCDF1_sumstats_munged.parquet \
    --non-funct \
    --n 2429 \
    --chr 15 \
    --start 26426847 \
    --end 26643821 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/CCDF1_NO_PRIORS.15.216kb.gz

python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.22 \
    --sumstats CCD_250822/geno_files/CCDF1_sumstats_munged.parquet \
    --non-funct \
    --n 2429 \
    --chr 22 \
    --start 25233990 \
    --end 25480179 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/CCDF1_NO_PRIORS.22.250kb.gz

python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.22 \
    --sumstats CCD_250822/geno_files/CCDF1_sumstats_munged.parquet \
    --non-funct \
    --n 2429 \
    --chr 22 \
    --start 25333990 \
    --end 25380179 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/CCDF1_NO_PRIORS.22.50kb.gz


#CCDF2
python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.1 \
    --sumstats CCD_250822/geno_files/CCDF2_sumstats_munged.parquet \
    --non-funct \
    --n 2506 \
    --chr 1 \
    --start 84225366 \
    --end 84429309 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/CCDF2_NO_PRIORS.1.204kb.gz

python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.1 \
    --sumstats CCD_250822/geno_files/CCDF2_sumstats_munged.parquet \
    --non-funct \
    --n 2506 \
    --chr 1 \
    --start 84325366 \
    --end 84329309 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/CCDF2_NO_PRIORS.1.4kb.gz

python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.5 \
    --sumstats CCD_250822/geno_files/CCDF2_sumstats_munged.parquet \
    --non-funct \
    --n 2506 \
    --chr 5 \
    --start 43413771 \
    --end 43652486 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/CCDF2_NO_PRIORS.5.240kb.gz


python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.5 \
    --sumstats CCD_250822/geno_files/CCDF2_sumstats_munged.parquet \
    --non-funct \
    --n 2506 \
    --chr 5 \
    --start 43513771 \
    --end 43552486 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/CCDF2_NO_PRIORS.5.40kb.gz



#CCDF3
python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.23 \
    --sumstats CCD_250822/geno_files/CCDF3_sumstats_munged.parquet \
    --non-funct \
    --n 2428 \
    --chr 23 \
    --start 27987645 \
    --end 28208643 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/CCDF3_NO_PRIORS.23.220kb.gz


#item151
python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.3 \
    --sumstats CCD_250822/geno_files/ITEM151_sumstats_munged.parquet \
    --non-funct \
    --n 2421 \
    --chr 3 \
    --start 55540298 \
    --end 55743099 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM151_NO_PRIORS.3.202kb.gz

python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.3 \
    --sumstats CCD_250822/geno_files/ITEM151_sumstats_munged.parquet \
    --non-funct \
    --n 2421 \
    --chr 3 \
    --start 55640298 \
    --end 55643099 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM151_NO_PRIORS.3.2kb.gz

#item152
python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.37 \
    --sumstats CCD_250822/geno_files/ITEM152_sumstats_munged.parquet \
    --non-funct \
    --n 2415 \
    --chr 37 \
    --start 29056938 \
    --end 29280887 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM152_NO_PRIORS.37.224kb.gz

#item150
python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.34 \
    --sumstats CCD_250822/geno_files/ITEM150_sumstats_munged.parquet \
    --non-funct \
    --n 2417 \
    --chr 34 \
    --start 6639027 \
    --end 6855846 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM150_NO_PRIORS.34.216kb.gz

python3 finemapper.py \
    --geno CCD_250822/geno_files/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.5 \
    --sumstats CCD_250822/geno_files/ITEM150_sumstats_munged.parquet \
    --non-funct \
    --n 2417 \
    --chr 5 \
    --start 88161974 \
    --end 88363559 \
    --method susie \
    --allow-missing \
    --max-num-causal 5 \
    --out CCD_250822/output/ITEM150_NO_PRIORS.5.202kb.gz