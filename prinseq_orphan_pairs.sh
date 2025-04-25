#!/bin/bash

# Define directories
INPUT_DIR="."             # Directory containing FastQ files
OUTPUT_GOOD="./orphan_pairs"  # Output folder for good files
OUTPUT_BAD="./orphan_pairs"   # Output folder for bad files (same as good files)
PRINSEQ_SCRIPT="./prinseq-lite-0.20.4/prinseq-lite.pl" # Path to Prinseq script

# Export the necessary variables for parallel to access
export PRINSEQ_SCRIPT
export OUTPUT_GOOD
export OUTPUT_BAD

# Use GNU Parallel to process paired FastQ files
ls ${INPUT_DIR}/*_filtered_1_cleaned.fastq | parallel -j 32 '
  FASTQ1={};
  FASTQ2="${FASTQ1/_filtered_1_cleaned.fastq/_filtered_2_cleaned.fastq}";

  # Run Prinseq to remove orphan pairs
  perl $PRINSEQ_SCRIPT \
    -fastq $FASTQ1 \
    -fastq2 $FASTQ2 \
    -min_len 30 \
    -out_good "${OUTPUT_GOOD}/$(basename $FASTQ1 _filtered_1_cleaned.fastq)_filtered_cleaned_new" \
    -out_bad "${OUTPUT_BAD}/$(basename $FASTQ1 _filtered_1_cleaned.fastq)_bad_new" \
    -out_format 3
'

