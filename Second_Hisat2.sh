#!/bin/bash

# Directories
sam_dir="/scratch/alice/a/ac876/HISAT2/Alignment_Output"
bam_dir="/scratch/alice/a/ac876/HISAT2/Alignment_Output/bam"
counts_dir="/scratch/alice/a/ac876/HISAT2/Alignment_Output/counts"
gtf_file="/scratch/alice/a/ac876/HISAT2/index/Homo_sapiens.GRCh37.87.gtf"
cpus=24

# Make sure output folders exist
mkdir -p "$bam_dir"
mkdir -p "$counts_dir"

# Process each SAM file
for sam_file in "$sam_dir"/*_filtered.sam; do
    # Extract the sample name before _filtered.sam
    sample_name=$(basename "$sam_file" "_filtered.sam")
    
    echo "Processing sample: $sample_name"

    # Convert SAM to BAM
    samtools view -bS "$sam_file" > "$bam_dir/${sample_name}.bam"

    # Sort BAM
    samtools sort "$bam_dir/${sample_name}.bam" -o "$bam_dir/${sample_name}.sorted.bam"

    # Index sorted BAM
    samtools index "$bam_dir/${sample_name}.sorted.bam"

    # Remove intermediate unsorted BAM to save space
    rm "$bam_dir/${sample_name}.bam"

    echo "Finished sample: $sample_name"
done

# Run featureCounts to create gene counts matrix
echo "Running featureCounts on all sorted BAM files..."

featureCounts -a "$gtf_file" -o "$counts_dir/gene_counts.txt" -t exon -g gene_id -p -T $cpus "$bam_dir"/*.sorted.bam

echo "Gene counts matrix created at $counts_dir/gene_counts.txt"

