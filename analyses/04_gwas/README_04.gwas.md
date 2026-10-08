# 04. GWAS
Running two models of gwas.

## 1_mlma-loco (factors):

Follow process for factor GWAS here: 

``` 
MLMA_PER_FACTOR.sh
``` 
Run clump in plink to define regions to finemap in 3_finemap_Susie: 

``` 
clump_plink_factors.sh
``` 

Export sumstats to use in 3_finemap_Susie:

``` 
export_gwas_sumstats_MLMA.R
``` 


## 2_polmm (items):

Follow process for item GWAS here. Run R script and change for each item:

``` 
POLMMgrab_PER_ITEM.sh
``` 
Run clump in plink to define regions to finemap in 3_finemap_Susie: 

``` 
clump_plink_items.sh
``` 

Export sumstats to use in 3_finemap_Susie:

``` 
export_gwas_sumstats_POLMM.R
``` 

## 3_finemap_Susie:

We used Susie without functional annotations to finemap clumped CCD associated loci, in the environment from polyfun: https://github.com/omerwe/polyfun


* ### To run SuSie finemap we need:

* #### 1. Make geno input: 
``` 
# one to use for all phenotypes, split by chr (e.g., CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.$chr)
finemap_create_genofiles.sh
```


* #### 2. Make sumstats harmonized and in parquet format, e.g., ITEM149_sumstats_munged.parquet

``` 
## 2A: harmonize sumstats: run py script for all phenotypes using this sbatch script:
2a_run_harmonize_sumstats_CCD.sbatch #uses harmonize_sumstats_CCD_SIZE.py OBS! change to your path

## 2B: munge sumstats into parquet format:
2b_munge_sumstats_parquet.sh #use polyfun py script to change the format of sumstats

```
* ### The scripts to run finemap on all CCD regions are here:
``` 
finemapNOFUNCT_multiple_CCDregions.sh 
finemapNOFUNCT_multiple_CCDregions2.sh 
finemapNOFUNCT_multiple_CCDregions3.sh 
finemapNOFUNCT_multiple_CCDregions4.sh 

##output is summarized in: (listing all SNPs in a credible set for each locus)
NO_PRIORS_finemap_SNPs_in_CS_251122.xlsx
``` 


## /data

* **geno files** - deposited genotype files are here... LINK!!!!
* **sumstat files** - deposited sumstat files are here... LINK!!!!

* **ALLFAM_DA_MERGED_GENCOVE_AXIOM_QC6.fam** - fam with dogs in all dogCD gwas:s
* **SIZE_DA_MERGED_GENCOVE_AXIOM_QC6.fam** - fam with dogs in SIZE gwas (Darwins Ark phenotype), part of finemapping pipeline in polyfun
* **STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.fam** - fam with dogs in STUCK gwas:s (Darwins Ark phenotype), part of finemapping pipeline in polyfun
* **NO_PRIORS_finemap_SNPs_in_CS_251122.xlsx** - output from finemapped dogCD loci, all SNPs in a Credible set.

