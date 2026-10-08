#!/bin/bash -l
#SBATCH -A $proj_id
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 35:00:00
#SBATCH -J harmonize_sumstats


source /sw/apps/conda/latest/rackham_stage/etc/profile.d/conda.sh
conda activate /home/kteng/.conda/envs/polyfun 

cd /path/to/

#python3 ../scripts/harmonize_sumstats_CCD_SIZE.py SIZE_sumstats.txt SIZE_sumstats_harmonized.txt
#python3 ../scripts/harmonize_sumstats_CCD_SIZE.py STUCK_sumstats.txt STUCK_sumstats_harmonized.txt
python3 ../scripts/harmonize_sumstats_CCD_SIZE.py CCDF1_sumstats.txt CCDF1_sumstats_harmonized.txt
python3 ../scripts/harmonize_sumstats_CCD_SIZE.py CCDF2_sumstats.txt CCDF2_sumstats_harmonized.txt
python3 ../scripts/harmonize_sumstats_CCD_SIZE.py CCDF3_sumstats.txt CCDF3_sumstats_harmonized.txt
python3 ../scripts/harmonize_sumstats_CCD_SIZE.py ITEM145_sumstats.txt ITEM145_sumstats_harmonized.txt
python3 ../scripts/harmonize_sumstats_CCD_SIZE.py ITEM146_sumstats.txt ITEM146_sumstats_harmonized.txt
python3 ../scripts/harmonize_sumstats_CCD_SIZE.py ITEM147_sumstats.txt ITEM147_sumstats_harmonized.txt
python3 ../scripts/harmonize_sumstats_CCD_SIZE.py ITEM148_sumstats.txt ITEM148_sumstats_harmonized.txt
python3 ../scripts/harmonize_sumstats_CCD_SIZE.py ITEM149_sumstats.txt ITEM149_sumstats_harmonized.txt
python3 ../scripts/harmonize_sumstats_CCD_SIZE.py ITEM150_sumstats.txt ITEM150_sumstats_harmonized.txt
python3 ../scripts/harmonize_sumstats_CCD_SIZE.py ITEM151_sumstats.txt ITEM151_sumstats_harmonized.txt
python3 ../scripts/harmonize_sumstats_CCD_SIZE.py ITEM152_sumstats.txt ITEM152_sumstats_harmonized.txt
python3 ../scripts/harmonize_sumstats_CCD_SIZE.py ITEM153_sumstats.txt ITEM153_sumstats_harmonized.txt
python3 ../scripts/harmonize_sumstats_CCD_SIZE.py ITEM154_sumstats.txt ITEM154_sumstats_harmonized.txt
python3 ../scripts/harmonize_sumstats_CCD_SIZE.py ITEM155_sumstats.txt ITEM155_sumstats_harmonized.txt
python3 ../scripts/harmonize_sumstats_CCD_SIZE.py ITEM7_sumstats.txt ITEM7_sumstats_harmonized.txt
python3 ../scripts/harmonize_sumstats_CCD_SIZE.py ITEM93_sumstats.txt ITEM93_sumstats_harmonized.txt
python3 ../scripts/harmonize_sumstats_CCD_SIZE.py ITEM95_sumstats.txt ITEM95_sumstats_harmonized.txt