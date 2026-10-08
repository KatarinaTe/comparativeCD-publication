#!/bin/bash

# Output file that will contain the combined results
output_file="combined_knownsites.depth.txt"

# Clear the output file if it already exists
> "$output_file"

# Loop over all matching files, regardless of whether the sample ID starts with
# letters, numbers, or a mix of both
for file in *.sorted.merged.MarkDups.BQSR.bam.knownsites.depth.txt; do
  # Skip the loop if no files match the pattern
  [ -e "$file" ] || continue

  # Remove the fixed suffix from the filename to get the sample ID
  sample_id=$(basename "$file" | sed 's/\.sorted\.merged\.MarkDups\.BQSR\.bam\.knownsites\.depth\.txt$//')

  # Extract the 3rd column from rows 2 and 3, then join them with tabs
  values=$(awk 'NR==2 || NR==3 {print $3}' "$file" | paste -sd '\t' -)

  # Write sample ID and values to the output file
  printf '%s\t%s\n' "$sample_id" "$values" >> "$output_file"
done