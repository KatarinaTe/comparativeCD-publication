#!/bin/bash

# Location of UU dog bed files

UU_CRE=/proj/caninebrain_2025/snic2022-6-164/brain_analysis/results/chromatin_states/ChromHMM/states9_bed_files

# Active states to retain: 
# Promoter-like: 1,4,5,6
# Enhancer-like: 2,8
# Open chromatin: 3

# Create beds with only those states to simplify the intersect

# for f in *.bed; do
#    # Extract tissue name (everything before the first underscore)
#    tissue=$(echo "$f" | cut -d'_' -f1)
#    
#    # Process file with awk
#    awk -v t="$tissue" 'BEGIN { FS="\t"; OFS="\t" } 
#        $4 ~ /^(1|4|5|6|2|8|3)$/ { print $1, $2, $3, $4, t }' "$f" > "${tissue}_UU_CRE_active_elements.bed"
# done

# Now test for overlaps with our GWAS regions
# Count bp overlap rather than elements as there seems to be a lot of elements in our dog brain data
# of variable length

#!/bin/bash

GWAS="/proj/caninebrain_2025/CCD_brain_CRE_intersect/CRE_overlap/size_regions/SIZE_gwas_clump100kb_regions_for_cCREintersect_clean.txt"
OUTPUT="SIZE_brain_region_bp_overlaps.txt"

echo -e "Region\tOverlap_bp\tTotal_bp" > "$OUTPUT"

REGIONS=("ACG" "cerebellum" "frontal" "hypothalamus" "occipital" "striatum" "temporal" "thalamus")

echo "Processing regions for Base Pair overlaps..."

for REGION in "${REGIONS[@]}"; do
    
    INPUT_BED="${REGION}_UU_CRE_active_elements.bed"

    # 1. Calculate Total background base pairs (N)
    # sum+0 ensures that if a file is empty, it prints "0" instead of a blank space
    TOTAL_BP=$(awk '{sum += $3 - $2} END {print sum+0}' "$INPUT_BED")

    # 2. Intersect and calculate Overlapping base pairs (O)
    # - Standard intersect outputs the exact overlapping segments
    # - We sort and merge to prevent double-counting in case any GWAS clumps overlap each other
    # - Finally, awk sums the lengths of these exact overlapping footprints
    OVERLAP_BP=$(bedtools intersect -a "$INPUT_BED" -b "$GWAS" | \
                 bedtools sort -i - | \
                 bedtools merge -i - | \
                 awk '{sum += $3 - $2} END {print sum+0}')

    # 3. Append the results to the output file
    echo -e "${REGION}\t${OVERLAP_BP}\t${TOTAL_BP}" >> "$OUTPUT"
    
    echo "Finished $REGION: $OVERLAP_BP bp overlaps out of $TOTAL_BP total bp."

done

echo "--------------------------------------------------"
echo "Done! Results saved to $OUTPUT"
echo "--------------------------------------------------"
