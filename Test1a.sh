# Define directories
index_dir="/scratch/alice/r/rsh31/Gene_Index"
input_dir="/scratch/alice/r/rsh31/Filtered_files"
output_dir="/scratch/alice/r/rsh31/STAR_align"


# Loop through all *_1_cleaned_new_ready.fastq files in the input directory
for file1 in $input_dir/*_1_cleaned_new_ready.fastq; do
  # Check if corresponding _2 file exists
  file2="${file1/_1_cleaned_new_ready.fastq/_2_cleaned_new_ready.fastq}"
  
  if [ -f "$file2" ]; then
    # Generate a prefix for output based on the input file name
    prefix=$(basename "$file1" _1_cleaned_new_ready.fastq)

    # Run STAR for paired-end alignment
    STAR --runThreadN 48 \
         --genomeDir $index_dir \
         --readFilesIn $file1 $file2 \
         --outFileNamePrefix $output_dir/${prefix}_ \
         --outSAMtype BAM SortedByCoordinate \
         --outSAMstrandField intronMotif \
         --alignSJoverhangMin 8 \
         --alignSJDBoverhangMin 1 \
         --outFilterMismatchNmax 999 \
         --outFilterMismatchNoverLmax 0.04 \
         --alignIntronMin 20 \
         --alignIntronMax 1000000 \
         --alignMatesGapMax 1000000 \
         --outFilterType BySJout \
         --outFilterMultimapNmax 20
         --quantMode GeneCounts 
  else
    # If the corresponding file2 does not exist, print a warning
    echo "Warning: corresponding file for $file1 not found ($file2). Skipping alignment for this pair."
  fi
done

