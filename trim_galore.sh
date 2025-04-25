#!/bin/bash

# Define directories
INPUT_DIR="."  # Current directory contains fastq files
OUTPUT_DIR="./trim_galore"  

# Use GNU Parallel to process paired FastQ files
ls ${INPUT_DIR}/*_filtered_cleaned_new_1.fastq | parallel -j 32 '
  FASTQ1={};
  FASTQ2="${FASTQ1/_filtered_cleaned_new_1.fastq/_filtered_cleaned_new_2.fastq}";

  # Run Trim Galore in paired-end mode
  ~/TrimGalore-0.6.10/trim_galore --paired --cores 8 \
    --stringency 1 \
    --nextera \
    --output_dir "$OUTPUT_DIR" \
    --fastqc \
    "$FASTQ1" "$FASTQ2"

  # Rename outputs to match naming convention
  mv "$OUTPUT_DIR/$(basename "$FASTQ1" _filtered_cleaned_new_1.fastq)_filtered_cleaned_new_1_val_1.fq" \
     "$OUTPUT_DIR/$(basename "$FASTQ1" _filtered_cleaned_new_1.fastq)_filtered_cleaned_new_ready_1.fastq"

  mv "$OUTPUT_DIR/$(basename "$FASTQ2" _filtered_cleaned_new_2.fastq)_filtered_cleaned_new_2_val_2.fq" \
     "$OUTPUT_DIR/$(basename "$FASTQ2" _filtered_cleaned_new_2.fastq)_filtered_cleaned_new_ready_2.fastq"
'

