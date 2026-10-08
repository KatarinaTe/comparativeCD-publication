#!/bin/bash -l
#SBATCH -A $proj_id
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 40:00:00
#SBATCH -J loco_mlmagrm_h2_CCDF1
#SBATCH --mail-user $user_id
#SBATCH --mail-type=ALL

module load bioinfo-tools gcta/1.94.1

cd /path/to/

#OBS! example is for Factor 1, do the same for factor 2-3:
#make grm
gcta64 --bfile CCDF1_DA_MERGED_GENCOVE_AXIOM_QC6 --autosome  --make-grm  --out CCDF1_DA_MERGED_GENCOVE_AXIOM_QC6 --autosome-num 38

#heritability
  echo "Performing GCTA-GREML with constraint..."
  gcta64 --grm CCDF1_DA_MERGED_GENCOVE_AXIOM_QC6 \
          --reml \
          --pheno CCDF1_DA_MERGED_GENCOVE_AXIOM_QC6.phen \
          --mpheno 1 \
          --covar sexCCDF1_DA_MERGED_GENCOVE_AXIOM_QC6.covar \
          --qcovar ageCCDF1_DA_MERGED_GENCOVE_AXIOM_QC6.qcovar \
          --thread-num 16 \
          --out CCDF1_QC6.REML.no-lds


  echo "Performing GCTA-GREML without constraint..."
    gcta64 --grm CCDF1_DA_MERGED_GENCOVE_AXIOM_QC6 \
            --reml \
            --reml-no-constrain \
            --pheno CCDF1_DA_MERGED_GENCOVE_AXIOM_QC6.phen \
            --mpheno 1 \
            --covar sexCCDF1_DA_MERGED_GENCOVE_AXIOM_QC6.covar \
            --qcovar ageCCDF1_DA_MERGED_GENCOVE_AXIOM_QC6.qcovar \
            --thread-num 16 \
            --out CCDF1_QC6.REML.no-lds.no-constraint 
          

#run gwas - loco:
gcta64 --mlma-loco --bfile CCDF1_DA_MERGED_GENCOVE_AXIOM_QC6 --pheno CCDF1_DA_MERGED_GENCOVE_AXIOM_QC6.phen --out CCDF1_QC6_LOCO --thread-num 16  --qcovar ageCCDF1_DA_MERGED_GENCOVE_AXIOM_QC6.qcovar --covar sexCCDF1_DA_MERGED_GENCOVE_AXIOM_QC6.covar  --autosome-num 38
#run gwas - with grm:
#gcta64 --mlma --bfile CCDF1_DA_MERGED_GENCOVE_AXIOM_QC6 --grm CCDF1_DA_MERGED_GENCOVE_AXIOM_QC6  --pheno CCDF1_DA_MERGED_GENCOVE_AXIOM_QC6.phen --thread-num 16 --qcovar ageCCDF1_DA_MERGED_GENCOVE_AXIOM_QC5.covar --covar sexCCDF1_DA_MERGED_GENCOVE_AXIOM_QC5.covar  --out CCDF1_QC6

Source  Variance        SE
V(G)    0.254735        0.071729
V(e)    0.777515        0.073378
Vp      1.032249        0.029674
V(G)/Vp 0.246776        0.068778
logL    -1269.094
logL0   -1282.918
LRT     27.647
df      1
Pval    7.2787e-08
n       2429
