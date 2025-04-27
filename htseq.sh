#!/bin/bash

# Define directories
bam_dir="/scratch/alice/s/sj485/srr_files/star_files/STAR_align"
gtf_file="/scratch/alice/s/sj485/srr_files/Homo_sapiens.GRCh37.87.gtf"
output_dir="/scratch/alice/s/sj485/srr_files/HTSeq_counts"

# Make output dir if it doesn’t exist
mkdir -p $output_dir

# Loop through all STAR BAMs
for bam in $bam_dir/*_Aligned.sortedByCoord.out.bam; do
    sample=$(basename "$bam" _Aligned.sortedByCoord.out.bam)

    htseq-count \
        -f bam \
        -r pos \
        -s no \
        -m intersection-nonempty \
        "$bam" "$gtf_file" > "$output_dir/${sample}_counts.txt"
done

