welcome! ʕ •ᴥ•ʔ 

PREPROCESSING - 

prinseq1.sh - Initial quality control and filtering on raw fastq file pairs

run_cutadapt_new.sh - Removal of overrepresented sequnces

prinseq_orphan_pairs.sh - To remove any new orphan pairs generated by the trimming process

trim_galore.sh - To remove Nextera adapters 

OLD PIPELINE - 

Test_1a.sh - The script tells STAR where to find your genome index, input FASTQ files, and output folder

htseq.sh - to quantify genes into counts matrix following alignment with STAR

gene_count_matrix.R - converting all gene count expression csvs into one gene count matrix

reanalysis.R - attempt to recreate original pipeline analysis

NEW PIPELINE - 

matrix.R- Converts original gene count CSVs into a Seurat-compatible gene matrix, to support familiarisation with the analysis workflow before working with the new gene matrix

Initial_Hisat2.sh - Initial alignment using alternative method (Hisat2) to get SAM files

Second_Hisat2.sh - Alignment and quantification using alternative method (Hisat2 and feature counts) to create new gene counts matrix 

AP_NEW1.R - Conversion and data analysis of new gene counts matrix using alternative pipeline (Seurat)


