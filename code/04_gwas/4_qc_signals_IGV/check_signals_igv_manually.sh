### We have bams from 1200 Darwin samples:

ll -h /proj/naiss2024-6-308/nobackup/private/Jennifer/DarwinsArkMapping/2C_mapping_Darwins1200/BamBaiOut > Darwin1200_bambai_sampleID.txt

awk '{split($NF, a, "."); print a[1]}' Darwin1200_bambai_sampleID.txt | sort | uniq > Darwin1200_bambai_sampleID_uniq.txt
awk 'NR > 1 {print $1, $1}' Darwin1200_bambai_sampleID_uniq.txt > Darwin1200_bambai_sampleID_uniq2.txt


### EXAMPLE REGION: check signal for Item149 chr1:117,641,155
/proj/snic2022-6-229/private/Darwin_Affy/NEW_DA_IMP_MERGED/2024-10-11/

grep chr1:117641155 geno_files/ITEM149_sumstats_harmonized.txt
#SNP     A1      A2      freq    BETA    se      P       N       CHR     POS
#chr1:117641155  A       G       0.0101156069364162      1.39512171847861        0.270706114445185       1.13382970960681e-07    4843.99999999999        1       117641155

plink --dog --bfile ../2024-10-11/DA_MERGED_GENCOVE_AXIOM_QC4 --snp chr1:117641155 --recode --out ../2024-10-11/leadSNP_item149chr1
grep A ../2024-10-11/leadSNP_item149chr1.ped


### I think I can just extract each SNP directly, it is quite obvious which is the rare allele:
#run for all 1200 ind for which we have bam-files for the lead SNPs of all associated loci:
#List SNP includes ALL snps in associated lead SNP list (including the ones with no genes or not part of a clump)
#SNP_selected.txt includes only lead SNPs contributing with a gene.

plink --dog --bfile ../../2024-10-11/DA_MERGED_GENCOVE_AXIOM_QC4 --keep Darwin1200_bambai_sampleID_uniq2.txt --extract SNP_selected.txt --recode --out Darwin1200_bambai_sampleID_leadSNP

#also check the regions not contributing with a gene to netcoloc:
plink --dog --bfile ../../2024-10-11/DA_MERGED_GENCOVE_AXIOM_QC4 --keep Darwin1200_bambai_sampleID_uniq2.txt --extract ../SNP_selected2.txt --recode --out Darwin1200_bambai_sampleID_leadSNP_secondround


# Working with IGV in Rackham - Sept 01 - 2025

1. Log in with -X
e.g., ssh -X <username>@rackham.uppmax.uu.se

2. Navigate to where the genome fasta is stored
/crex/proj/uppstore2017228/KLT.03.CGEN/Mischka/Other_dog_genome_assembly/GSD1.0_RosY/UU_Cfam_GSD_1.0_ROSY.fa

3a. To avoid black screen run this before launching igv-node:
export _JAVA_OPTIONS='-Dsun.java2d.xrender=false -Dsun.java2d.pmoffscreen=false'

3. Load IGV module
module load bioinfo-tools IGV
igv-node

4. In the IGV GUI, make the CF4 genome
Genomes -> Load Genome from file -> 
Choose
/crex/proj/uppstore2017228/KLT.03.CGEN/Mischka/Other_dog_genome_assembly/GSD1.0_RosY/UU_Cfam_GSD_1.0_ROSY.fa

5. Upload the bams you want to look at. [/proj/snic2022-6-229]
File -> Load from file ->
/proj/naiss2024-6-308/nobackup/private/Jennifer/DarwinsArkMapping/2C_mapping_Darwins1200/BamBaiOut

You may have to look for [/proj/snic2022-6-229]