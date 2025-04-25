#!/bin/bash

# Directory where the fastq files are located
fastq_dir="/scratch/alice/a/ac876/new_files/test_output_Prinseq/orphan_pairs/trim_galore"
# HISAT2 index directory
index_dir="/scratch/alice/a/ac876/HISAT2/index"
# Output directory 
output_dir="/scratch/alice/a/ac876/HISAT2/Alignment_Output"
# Annotation GTF file 
gtf_file="/scratch/alice/a/ac876/HISAT2/index/Homo_sapiens.GRCh37.87.gtf"

# Create the output directory 
mkdir -p $output_dir

# Number of CPUs to use for parallelisation
CPUs=24

# Loop through each pair of fastq files (they follow the pattern *_1_cleaned_new_ready.fastq and *_2_cleaned_new_ready.fastq)
for R1 in $fastq_dir/*_1_cleaned_new_ready.fastq; do
    # Get the R2 file by replacing the _1 part with _2
    R2="${R1/_1_cleaned_new_ready.fastq/_2_cleaned_new_ready.fastq}"
    
    # Extract sample name from the R1 filename (before the _1_cleaned_new_ready.fastq part)
    sample_name=$(basename "$R1" "_1_cleaned_new_ready.fastq")

    echo "Processing sample: $sample_name"

    # Align the paired-end reads with HISAT2, using parallelisation (-p)
    hisat2 -x /scratch/alice/a/ac876/HISAT2/index/hg19_index/genome -1 $R1 -2 $R2 -S $output_dir/"$sample_name".sam -p $CPUs


    echo "Finished processing sample: $sample_name"
done

echo "All samples processed!"
