#!/bin/bash

# Loop over all _filtered_1.fastq files in the directory and run Cutadapt in parallel
find . -name "*filtered_1.fastq" | parallel -j 32 --bar '
    # Get the base name (remove the "_filtered_1.fastq" part)
    base_name=$(echo {} | sed "s/_filtered_1.fastq//")

    # Run Cutadapt for each pair of files
    cutadapt \
        -a file:overrepresented_sequences.fasta \
        -b file:overrepresented_sequences.fasta \
        -m 30 \
        -e 0.15 \
        -o "${base_name}_filtered_1_cleaned.fastq" \
        -p "${base_name}_filtered_2_cleaned.fastq" \
        "${base_name}_filtered_1.fastq" \
        "${base_name}_filtered_2.fastq"
'

