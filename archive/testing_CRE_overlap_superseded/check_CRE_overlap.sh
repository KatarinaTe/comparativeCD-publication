#!/bin/bash

# Location of EPIC dog bed files

EPIC=/proj/canine_compgen/uppstore2017228/KLT.04.200M/200m_MC/dog_brain_annotation/EPIC_dog/chromatin_states/CanFam4

#Epic states:

#state  type
#1       Bivalent TSS/Enh
#2       Repressed polycomb
#3       Repressed
#4       Heterochromatin
#5       ZNF genes + repeats
#6       Weak TSS
#7       Flanking active TSS1
#8       Flanking active TSS2
#9       Active TSS
#10      Quiescent
#11      Active, weak enhancer
#12      Active, strong enhancer
#13      Active, poised enhancer

#So we want to intersect with states 7-9 (promoter-like) and 11-13 (enhancer-like)

# Create beds with only those states to simplify the intersect

#for f in *.bed; do
#    # Extract tissue name (everything before the first underscore)
#    tissue=$(echo "$f" | cut -d'_' -f1)
#    
#    # Process file with awk
#    awk -v t="$tissue" 'BEGIN { FS="\t"; OFS="\t" } 
#        $4 ~ /^(7|8|9|11|12|13)$/ { print $1, $2, $3, $4, t }' "$f" > "${tissue}_EPIC_active_elements.bed"
#done

# Pool and collapse brain regions into a single file

# Define the GWAS file
GWAS="gwas_clump100kb_regions_for_cCREintersect.txt"

echo "Step 1: Pooling, sorting, and merging Brain tissues..."
# Concatenate cerebellum and cerebrum, sort by coordinates, and merge overlapping regions
cat cerebellum_EPIC_active_elements.bed cerebrum_EPIC_active_elements.bed | \
    bedtools sort -i - | \
    bedtools merge -i - > EPIC_active_merged_brain.bed

echo "Step 2: Pooling, sorting, and merging Other tissues..."
# Concatenate the remaining 9 tissues
cat colon_EPIC_active_elements.bed kidney_EPIC_active_elements.bed liver_EPIC_active_elements.bed \
    lung_EPIC_active_elements.bed mammary_gland_EPIC_active_elements.bed ovary_EPIC_active_elements.bed \
    pancreas_EPIC_active_elements.bed spleen_EPIC_active_elements.bed stomach_EPIC_active_elements.bed | \
    bedtools sort -i - | \
    bedtools merge -i - > EPIC_active_merged_other.bed

echo "Step 3: Calculating Background Totals (N)..."
# Count the total number of distinct elements in each category
N_brain=$(wc -l < EPIC_active_merged_brain.bed)
N_other=$(wc -l < EPIC_active_merged_other.bed)

echo "Step 4: Intersecting with GWAS regions (O)..."
# Use -u to report each tissue element only once, even if it overlaps multiple GWAS clumps
bedtools intersect -a EPIC_active_merged_brain.bed -b "$GWAS" -u > overlap_brain.bed
bedtools intersect -a EPIC_active_merged_other.bed -b "$GWAS" -u > overlap_other.bed

# Count the number of overlapping elements
O_brain=$(wc -l < overlap_brain.bed)
O_other=$(wc -l < overlap_other.bed)

echo "--------------------------------------------------"
echo "Done! Here is your Fisher's Exact Test contingency table:"
echo "--------------------------------------------------"
echo -e "\t\tOverlap with GWAS\tNo Overlap with GWAS"
echo -e "Brain\t\t$O_brain\t\t\t$((N_brain - O_brain))"
echo -e "Other\t\t$O_other\t\t\t$((N_other - O_other))"
echo "--------------------------------------------------"

######
#
# Ok, now remove the +-100kb 'tails' and test again on just the clump regions

#!/bin/bash

# INPUT_BED="gwas_clump100kb_regions_for_cCREintersect.txt"
# OUTPUT_BED="gwas_clump_original.bed"
# 
# echo "Shrinking regions by 100kb on both sides..."
# 
# awk -v shrink=100000 'BEGIN { OFS="\t" } {
#     # Add 100kb to the start coordinate
#     $2 = $2 + shrink
#     
#     # Subtract 100kb from the end coordinate
#     $3 = $3 - shrink
#     
#     # Sanity check: ensure the region is still valid (start is before end)
#     if ($2 < $3) {
#         print $0
#     } else {
#         # Print a warning to the console if a region breaks
#         print "Warning: Region on line " NR " (" $1 ") became invalid (start >= end) and was skipped." > "/dev/stderr"
#     }
# }' "$INPUT_BED" > "$OUTPUT_BED"
# 
# echo "Done! Shrunk regions saved to $OUTPUT_BED"
