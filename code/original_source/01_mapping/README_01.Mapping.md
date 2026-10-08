# 01. Mapping
Handles the initial processing of sequencing data for both high and low-coverage samples.

## 1_High-LowPass Mapping: 

Following instructions from: https://github.com/jmkidd/dogmap and see mapping script origin: https://github.com/Chao912/dog_10k/blob/main/dog_10k_mapping.sh

#### First, within the dog10k_mapping.sh script update links to the reference files:
``` 
#all of these ref genome files are found here: https://github.com/jmkidd/dogmap
assembly_ref=/path/to/UU_Cfam_GSD_1.0_ROSY.fa
known_variants=/path/to/UU_Cfam_GSD_1.0.BQSR.DB.bed.gz
depth_sites=/path/to/SRZ189891_722g.simp.header.CanineHD.names.GSD_1.0.filter.vcf.gz

#additional paths needed (see https://github.com/Chao912/dog_10k/blob/main/dog_10k_mapping.sh)
chunks_path=/path/to/chunks_10/ #split the genome into 10 chunks  
bwa2_path=/path/to/bwa-mem2-2.1_x64-linux/
gatk_path=/path/to/gatk-4.2.0.0/

#available in the 1_High.LowPass_mapping folder:
run_stats=/path/to/run-stats.py
``` 

### HighPass:
#### 1. Run mapping script:

``` 
#Use your sample list (in /01_mapping/data) and split:
split -l 1 High_pass_samples.txt SRR 

#Submit the job with script
for i in SRR*;do echo $i; sbatch dog10k_mapping.sh $i; done
``` 
#### 2. Prepare cohort.sample_map, see file ./data folder, to use in 1_HighPass_SNPfiltering:

``` 
SRR13340513 /path/to/SRR13340513/SRR13340513.sorted.merged.MarkDups.BQSR.g.vcf.gz
etc.
``` 

### LowPass: 

* ***OBS! 12 samples were included in these processes but excluded in the end (02_data_processing/3_Merging_filtering/01_gencove_axiom_merging_filtering.sh) due to no Darwin's Ark dog ID and no phenotypes. These are not available in SRA. These samples are listed here: 01_mapping/data/exclude_12_LowPass_samples.txt***

#### 1. Run mapping script:

``` 
#Use your sample list (in /01_mapping/data) and split:
split -l 1 All.DA_LowPass.Samples_fastq.txt DA 

#Submit the job with script
for i in DA*;do echo $i; sbatch dog10k_mapping.sh $i; done
``` 

#### 2. Extract depth from all samples from the stats run within the dog10k_mapping.sh:

``` 
extract_depth_autosomes_and_x.sh
``` 


#### 2. Move bam output to a .BamBaiOut folder:

``` 
./bamMove.sh All.DA_LowPass.Samples_Movelist
``` 

## 2_HighPass_SNPfiltering: 
Joint calling and SNP hard filtering (`vcf` output). And merging chr files in the end. Resulting in  a SNP.HF.ann.id.vcf.gz, SNP.HF.ann.id.vcf.gz.tbi 

####  run this with your cohort.sample_map in the same folder:
``` 
genomicsDB_joincall_SNPIDandHardFilter.sh
``` 

## 3_LowPass Imputation: 
Reference panel preparation and imputation using **GLIMPSE v1.1** to handle low-depth samples. The flow of the process can be found by following  https://odelaneau.github.io/GLIMPSE/glimpse1/tutorial_b38.html

Input = bam files

#### Follow process here:

``` 
LowPass_imputation.md
``` 

## /data

* **High_pass_samples.txt** - 11 high-pass samples, data available in SRA ""LINK"" OBS! update paths.
* **All.DA_LowPass.Samples_fastq.txt** - 3044 low-pass samples, data available in SRA ""LINK"" OBS! update paths.
* **All.DA_LowPass.Samples_bam.txt** - 3044 low-pass samples, links to bam files
* **All.DA_LowPass.Samples_MoveList** - 3044 low-pass sample IDs
* **CHR1-38.txt** - simple list of dog autosomes
* **cohort.sample_map** - HighPass sample list for input to SNP-filtering (sampleID, path/to/*g.vcf.gz)
* **exclude_12_LowPass_samples.txt** - 12 samples were included in the mapping processes but excluded in the end due to no Darwin's Ark dog ID and no phenotypes. These are not available in SRA. 
