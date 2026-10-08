#################################
## Low Pass Imputation from bam ##
#################################

### The flow of the process can be found by following 
https://odelaneau.github.io/GLIMPSE/glimpse1/tutorial_b38.html

### Running Imputation

**"Remember to update all paths & Sample naming"**

### Folders 
Suggestion is to create these folder setup during the process:
``` 
/0_Scripts
/1_phased-imputation-panel
/2_mapping
/3_GenotypeLikelihood
/4_ImputeAndPhase
/5_LigatedImputedCHRs 
/6_Final_imputed_dataset
```

### For GLIMPSE ->
```
path/to/GLIMPSEv1.1
```

### Scripts & Lists
```
01_VariableSites.sh
02_GenotypeLiklihood.sh
03_ImputeAndPhaseChrsByChunks.sh
04_LigatePerChromsome.sh
05_QCfilter_imputed.sh
06_vcf_to_plink.sh
07_Merge_chr_plink.sh

DoOnce_GlimpseChunks.sh
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

* **Download phased-impuation-panel here:**

https://kiddlabshare.med.umich.edu/dog10K/phased-imputation-panel/

``` 
AutoAndXPAR.Dog10K.Phased.bcf
AutoAndXPAR.Dog10K.Phased.bcf.csi
```

#### All samples are mapped with bams and bais here:
``` 
/2_mapping/BamBaiOut/ 
```

* ### Step 01 - extracting the polymorphic, biallelic SNP from the phased Dog10K panel. See
```
/0_Scripts/01_VariableSites.sh
```

* ### Step 02 - Calculate the GenotypeLikihoods of all samples by Chromosome -> THIS IS THE LONGEST STEP can take ~2 days for chr1 on 1800 samples
``` 
/0_Scripts/02_GenotypeLiklihood.sh
``` 
#### Check that paths are correct for all the samples and the script is run by 
```
./02_GenotypeLiklihood.sh CHR1-38.txt 
```

* ### Step 03 - Impute based on GenotypeLilihoods -> Around 7hours per chr. See 
```
03_ImputeAndPhaseChrsByChunks.sh
```

#### Here the chucks were pre-defined by the GLIMPSE software, and each chr is broken into many overlapping chunks
#### The input for the script is one line, with those chunks defined.
#### Generating the chunk input files you run for example CHR1
```
split -l 1 path/to/GLIMPSEv1.1/GLIMPSE_chr_chunks/chunks_chr1.txt chr1 
```

#### We group them in six folders (by running 10 chrs at a time, we allow for errors to be found):

* #### Make folders for each Chunk group under /0_Scripts: i.e. 
ChunkList.Chr1 

ChunkList.Chr2 

ChunkList.Chr3-10 

ChunkList.Chr11-20 

ChunkList.Chr21-30 

ChunkList.Chr31-37 

ChunkList.Chr38

* #### In each folder run Split.sh, e.g. for ChunkList3-10:

```
#!/bin/bash -l
#SBATCH -A $proj_id
#SBATCH -p core
#SBATCH -J Split
#SBATCH -t 10:00

module load bioinfo-tools
module load bcftools
module load samtools

cd /path/to/0_Scripts/ChunkList.Chr3-10

split -l 1 /path/to/GLIMPSEv1.1/GLIMPSE_chr_chunks/chunks_chr3.txt chr3
split -l 1 /path/to/GLIMPSEv1.1/GLIMPSE_chr_chunks/chunks_chr4.txt chr4
split -l 1 /path/to/GLIMPSEv1.1/GLIMPSE_chr_chunks/chunks_chr5.txt chr5
split -l 1 /path/to/GLIMPSEv1.1/GLIMPSE_chr_chunks/chunks_chr6.txt chr6
split -l 1 /path/to/GLIMPSEv1.1/GLIMPSE_chr_chunks/chunks_chr7.txt chr7
split -l 1 /path/to/GLIMPSEv1.1/GLIMPSE_chr_chunks/chunks_chr8.txt chr8
split -l 1 /path/to/GLIMPSEv1.1/GLIMPSE_chr_chunks/chunks_chr9.txt chr9
split -l 1 /path/to/GLIMPSEv1.1/GLIMPSE_chr_chunks/chunks_chr10.txt chr10
```

* #### In each folder, make the 03_ImputeAndPhaseChrsByChunks.sh adapted for each Chunk group and run the script like this, e.g. for ChunkList3-10: 
```
for i in chr*;do echo $i; sbatch 03_ImputeAndPhaseChrsByChunks_Chr3-10.sh $i; done
```

#### Note! for the script to run properly, your split line inputs need to match the path in line 58 of the script, i.e., done < $ScriptsFolder/$input_file

* ### Step 04 - Ligate the chunks you made above together. See per chr script -> 
```
04_LigatePerChromsome.sh
```

#### you will need to supply a list of the chunks to ligate in order. Keep them per chr to make the filtering on INFO field faster

* ### Step 05 - QC filter on imputed chr 
```
05_QCfilter_imputed_chr.sh
```
#### Note! you will need to remove the inter mediate bcf and vcf files. The final output is a chr*.merged.vcf.gz and chr*.merged.vcf.gz.csi 

* ### Step 06 - make vcf to plink files 
```
06_vcf_to_plink.sh
```

* ### Step 07 - merge chr 
```
07_merge_chr_plink.sh
```

