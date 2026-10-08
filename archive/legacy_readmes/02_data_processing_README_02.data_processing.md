# 02. Data processing
Integrates genotype array data with phenotypic surveys.

## 1_Survey Data:

Factor analysis/IRT for survey phenotypes. See:
``` 
Darwin_dog-CD_3-factorsolution.R
``` 
Includes plotting of histogram of the Factor 1-3 distributions for all survey dogs.

## 2_Axiom Imputation:

LiftOver operations (Cf3 to Cf4) and imputations. See scripts:

``` 
01_axiom_liftover.sh

02_axiom_impute.sh

03_axiomimp_to_plink.sh
``` 


#### Dog ref (UU_Cfam_GSD_1.0_ROSY) and phased imputation panel (dog10K)

* **Download reference files from here:**

https://zenodo.org/records/8084059

Jeffrey Kidd. (2023). Genome sequencing of 2,000 canids advances the understanding of demography, genome function and architecture [Data set]. Zenodo. https://doi.org/10.5281/zenodo.8084059

``` 
UU_Cfam_GSD_1.0_ROSY.fa
UU_Cfam_GSD_1.0.BQSR.DB.bed.gz
SRZ189891_722g.simp.header.CanineHD.names.GSD_1.0.filter.vcf.gz
``` 

* **Download phased-imputation-panel here:**

https://kiddlabshare.med.umich.edu/dog10K/phased-imputation-panel/

``` 
AutoAndXPAR.Dog10K.Phased.bcf
AutoAndXPAR.Dog10K.Phased.bcf.csi

#make one ref per chr: DO THIS ONCE #done within 02_axiom_impute.sh 
bcftools view /path/to/AutoAndXPAR.Dog10K.Phased.bcf -r $chrN  -Oz -o $REF_PATH'modi_AutoAndXPAR.Dog10K.Phased_'$chrN'.vcf.gz'
bcftools index $REF_PATH/'modi_AutoAndXPAR.Dog10K.Phased_'$chrN'.vcf.gz'

```

## 3_Merging & Filtering:

Merging the imputed Axiom and Gencove low-pass genotype data.
Filtering: Quality control for individuals (duplicates, phenos, and relatedness) and sample quality (depth >0.3) and genotypes (HWE, MAF).

``` 
01_gencove_axiom_merging_filtering.sh

02_create_ALLFAM_check_stats.sh #check all individuals left after QC for downstream analyses

03_create_files_per_factor_mlma.R

04_create_files_per_item_polmm.R
``` 
Includes processes back and forth between R and plink (terminal). Resulting in final plink files for each phenotype (N=17)

Includes plotting of histogram of the Factor 1-3 and items distributions for all genotyped dogs.

## 4_Population structure:

Population Structure: PCA generation using plink2.



## /data


1_Survey_data: 
* **response_df_dogCDitems.txt** - 1_Survey_data: includes the responses for 14 dogCD items from 25,302 dogs (from the Darwin's Ark survey), OBS! NA was omitted from EFA -> 12,989 dogs

2_Axiom_imputation:
* **ind_to_keep_round3.txt** - 2_Axiom_imputation: 411 axiom array samples
* **chr4_chr_start_stop.bed** - 2_Axiom_imputation: input files used in imputation
* **chromosomes.txt** - 2_Axiom_imputation: input files used in imputation

3_Merging_filtering:
* **DA_MERGED_GENCOVE_AXIOM_QC3modi.fam** 
* **DarwinsDogs_Q121_height_age_sex_batch_GENCOVE_AXIOM_QC4.tsv** 
* **DarwinsDogs_N-3285_bam_meandepths.txt** 
* **data6_2024-10-14.txt**  
* **Duplicates_to_keep_DarwinsDogs_N-9_21samples_bam_meandepths.txt** - 
* **fam_2584dogs.txt** - all dogs included in any of the gwas:s
* **list_exclude.txt** - Excluding duplicates and related dogs
* **list_exclude_secondexclusion.txt** - Excluding duplicates and related dogs
* **list_exclude_thirdexclusion.txt** - Excluding duplicates and related dogs
