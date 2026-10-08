#!/bin/bash
# full path to the outputfile
output_file=$1
# path to the samtools coverage output file
bamcov_dir=$2

# Clear the output file if it exists
> ${output_file}

cd ${bamcov_dir}

for file in [0-9]*; do

    # Extract sample ID from filename (this batch just has numbers as dog IDs in the file names)
    sample_id=$(echo "${file}" | grep -oP '^\d+')
    echo "${sample_id}"
     # Calculate mean of column 7. rows 1-38 (excluding header)
    mean=$(awk 'NR > 1 && NR <= 38 {sum += $7; count++} END {if (count > 0) print sum / count; else print "N/A"}' ${file})

    # Combine sample ID with mean. then append to output file
    echo -e "${sample_id}\t${mean}" >> ${output_file}

done
