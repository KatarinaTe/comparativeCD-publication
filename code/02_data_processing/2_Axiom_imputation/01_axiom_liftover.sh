########################################################################
########################################################################
#Axiom liftover cf3 to cf4 then imputation dog10K
########################################################################
########################################################################
#2024-10-23

##Data is here LINK.....!!!!
###ind_to_keep_round3.txt is in ./data

plink --bfile Affy_merged --dog --make-bed --keep ind_to_keep_round3.txt --out affy_round3
#1011992 variants and 411 dogs pass filters and QC.


####sex check using plink on AXIOM####
/Users/katst103/Documents/Dog-group/2020/Autism_Border_Collie/Morrill_data/canfam4_remapping/Axiom_remapping/affy_round3

#Note that, when there are n autosome pairs, the X chromosome is assigned numeric code n+1, Y is n+2, XY (pseudo-autosomal region of X) is n+3, and MT (mitochondria) is n+4.
#there is no X-chr....?

38	chr38:23902151	0	23902151	A	G
38	chr38:23903967	0	23903967	T	C
38	chr38:23910486	0	23910486	A	G
41	chr41:11685	0	11685	G	A
41	chr41:74057	0	74057	0	C

41 = XY (pseudo-autosomal region of X) is n+3
976991 - 1011992 #35001 markers

plink --dog --bfile affy_round3 --split-x 6600000 123798852 --indep-pairphase 20000 2000 0.5  --check-sex 0.5 0.8 --out sex_affy_round3



##########################################################################
#The liftover process done locally was as follows: NEW 2024-02-08
#########################################################################
#We used LiftOver command-line program (LINUX 64-bit x86, 2022-01-31)
#and canFam3ToCanFam4.over.chain.gz (2020-05-12) 

###use axiom plink files:
/path/to/Affy_merged.bed #- canfam3: is the axiom genotype data from Darwins Ark (array A+B) mapped to canfam3, 804 individuals
/path/to/affy_round3.fam ##- canfam3: consists of 411 individuals (specified in ind_to_keep_round3.txt)

#I work in this folder 
/path/to/axiom_remapping

#make correct format of chr and position etc.
awk '{print "chr"$1,$4-1,$4,$2,"0","+"}' affy_round3.map > affy_round3_cf3_BED.bed ## 

#run liftover
./liftOver affy_round3_cf3_BED.bed canFam3ToCanFam4.over.chain output.bed unlifted.bed

unlifted.bed #74202 -> 37101 variants not lifting
awk '{print $4}' output.bed > lifted_markers.txt #974891

#flip the ones coded as -
grep - output.bed > flip_markers.txt #47707
awk '{print $4}' flip_markers.txt > flip_markers_list.txt
   
   
#run plink in canfam3 and keep these markers
plink --bfile ../affy_round3 --extract lifted_markers.txt --flip flip_markers_list.txt --geno 0.05 --maf 0.01 --hwe 'midp' 0.000000000000001  --make-bed --dog --out affy_cf3 --snps-only 'just-acgt'
#--flip: 47670 SNPs flipped
#68377 variants removed due to missing genotype data (--geno).
#37 variants had at least one non-A/C/G/T allele name
#--hwe: 3738 variants removed due to Hardy-Weinberg exact test.
#156001 variants removed due to minor allele threshold(s)
#746059 variants and 411 dogs pass filters and QC.

#check if any variants do not have two specified alleles: NO
cat affy_cf3.bim | awk '{if ($5 == "0") print $0;}' > A1_is_0.txt #0
cat affy_cf3.bim | awk '{if ($6 == "0") print $0;}' > A2_is_0.txt #0
cat affy_cf3.bim | awk '{if ($5 == !"A/C/T/G") print $0;}' > A1_notACTG.txt #0
awk 'length($5) > 1 || length($6) > 1' affy_cf3.bim > indels.txt #0
 
####canfam4_build - SNP and pos ####

#exclude markers that lifted to another chromosome 
awk -F'\t' '{split($4, a, ":"); if ($1 != a[1]) print}' output.bed > no_match.txt #530
awk '{print $4}' no_match.txt > no_match_markers_remove.txt #530

plink --bfile affy_cf3 --exclude no_match_markers_remove.txt --dog --make-bed --out affy_cf3_a
#745854 variants and 411 dogs pass filters and QC.

#update map positions
awk '{print $4,$3}' output.bed > map_canfam4_build.txt #chr1	173603	173604	chr1:212740	0	+ --> it is the third column that has my new coordinate = 173604 (cf4) 212740 (cf3)
awk '{print $4,$1}' output.bed > chr_canfam4_build.txt
plink --bfile affy_cf3_a --update-map map_canfam4_build.txt --update-chr chr_canfam4_build.txt --dog --make-bed --out affy_cf4
#745854 variants and 411 dogs pass filters and QC.
#Warning: Base-pair positions are now unsorted!

#run again to sort the positions
plink --bfile affy_cf4 --dog --make-bed --out affy_cf4b

###update SNPID###
#update variant name: canfam4_build - Old_SNPID, newSNPID
awk '{print $2,"chr"$1":"$4}' affy_cf4b.bim > newSNPID.txt
plink --bfile affy_cf4b --update-name newSNPID.txt --dog --make-bed --out DA_AFFY
#Total genotyping rate is 0.995398.
#745854 variants and 411 dogs pass filters and QC.

