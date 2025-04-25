#!/bin/bash

# Define directories
INPUT_DIR="./fastq_out2"             # Directory containing FastQ files
OUTPUT_GOOD="./test_output_Prinseq" # Output folder for good files
OUTPUT_BAD="./test_output_Prinseq"  # Output folder for bad files (same as good files)
PRINSEQ_SCRIPT="./prinseq-lite-0.20.4/prinseq-lite.pl" # Path to Prinseq script

# Export the necessary variables for parallel to access
export PRINSEQ_SCRIPT
export OUTPUT_GOOD
export OUTPUT_BAD

# Use GNU Parallel to process paired FastQ files
ls ${INPUT_DIR}/*_1.fastq | parallel -j 32 '
  FASTQ1={};
  FASTQ2="${FASTQ1/_1.fastq/_2.fastq}";

  # Run Prinseq
  perl $PRINSEQ_SCRIPT \
    -fastq $FASTQ1 \
    -fastq2 $FASTQ2 \
    -min_len 30 \
    -trim_left 10 \
    -trim_qual_right 25 \
    -lc_method entropy \
    -lc_threshold 65 \
    -out_good "${OUTPUT_GOOD}/$(basename $FASTQ1 _1.fastq)_filtered" \
    -out_bad "${OUTPUT_BAD}/$(basename $FASTQ1 _1.fastq)_bad" \
    -out_format 3
'

