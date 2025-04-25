# List all .csv files in the directory https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/list.files
file_list <- list.files(pattern = "\\.csv$") 

# Read each file into a list of data frames https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/lapply
# Useful youtube video: https://www.youtube.com/watch?v=34sbvhr_pm8&t=145s&pp=ygUaZGF0YWRhZnQgaG93IHRvIHVzZSBsYXBwbHk%3D
gene_counts_list <- lapply(file_list, function(file) { 
  read.table(file, header = TRUE, row.names = 1) 
}) 

# Combine the list of data frames into one matrix (columns = cells) https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/do.call
gene_counts_matrix <- do.call(cbind, gene_counts_list) 

# Convert to numeric matrix (important for Seurat) https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/matrix
gene_counts_matrix <- as.matrix(gene_counts_matrix) 

# Handle duplicate gene names https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/duplicated
# Make.unique(): https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/make.unique
if (any(duplicated(rownames(gene_counts_matrix)))) {
  rownames(gene_counts_matrix) <- make.unique(rownames(gene_counts_matrix))
}

# Extract cell names from file names https://www.digitalocean.com/community/tutorials/sub-and-gsub-function-r
colnames(gene_counts_matrix) <- sub("\\.csv$", "", file_list)


