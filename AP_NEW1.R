# Load libraries
library(Seurat)  # Single-cell RNA-seq analysis
library(ggplot2)  # Data visualization
library(dplyr)  # Data manipulation
library(biomaRt)  # Biological data queries
library(patchwork)  # Combine ggplot2 plots
library(presto)  # Fast differential expression & clustering


# Gene counts matrix #

# Read the single gene matrix CSV (no header, tab-separated)
gene_counts_matrix <- read.table("gene_counts.csv", header = FALSE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)

# Assign default column names e.g. Sample1, Sample2, ...
colnames(gene_counts_matrix) <- paste0("Sample", seq_len(ncol(gene_counts_matrix)))

# To ensure it's a numeric matrix
gene_counts_matrix <- as.matrix(gene_counts_matrix)

# Handle duplicate gene names (Ensembl) just in case
if (any(duplicated(rownames(gene_counts_matrix)))) {
  rownames(gene_counts_matrix) <- make.unique(rownames(gene_counts_matrix))
}


# Create Seurat object with Ensembl IDs (original)
seurat_object_ensembl <- CreateSeuratObject(counts = gene_counts_matrix)
seurat_object_ensembl <- NormalizeData(seurat_object_ensembl, normalization.method = "LogNormalize", scale.factor = 1e6)


# Map Ensembl IDs to Gene Symbols using biomaRt
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_annotations <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = rownames(gene_counts_matrix),
  mart = ensembl
)

# Filter out rows where HGNC symbol is blank
gene_annotations <- gene_annotations[gene_annotations$hgnc_symbol != "", ]
# Match Ensembl IDs to your matrix
matched <- gene_annotations[match(rownames(gene_counts_matrix), gene_annotations$ensembl_gene_id), ]
# Replace Ensembl IDs with HGNC symbols
rownames(gene_counts_matrix) <- make.unique(matched$hgnc_symbol)
# Remove unmapped rows
gene_counts_matrix <- gene_counts_matrix[!is.na(rownames(gene_counts_matrix)) & rownames(gene_counts_matrix) != "", ]


# Create Seurat object with HGNC gene symbols
seurat_object_symbols <- CreateSeuratObject(counts = gene_counts_matrix)
seurat_object_symbols <- NormalizeData(seurat_object_symbols, normalization.method = "LogNormalize", scale.factor = 1e6)



#  Add Metadata from .soft File #

# Read in the .soft file as lines of text
soft_file <- readLines("GSE67835_family.soft")

# Split the file into blocks by individual sample entries
sample_blocks <- split(soft_file, cumsum(grepl("^\\^SAMPLE", soft_file)))

# Extract metadata (Sample ID, developmental stage, and cell type) from each sample block
sample_metadata <- lapply(sample_blocks, function(block) {
  id_line <- grep("^\\^SAMPLE = ", block, value = TRUE)  # Extract line containing GSM ID
  stage_line <- grep("!Sample_characteristics_ch1 = age:", block, value = TRUE)  # Extract developmental stage
  type_line <- grep("!Sample_characteristics_ch1 = cell type:", block, value = TRUE)  # Extract cell type
  
  # Clean the extracted values
  gsm <- if (length(id_line) > 0) sub("^\\^SAMPLE = ", "", id_line) else NA
  stage <- if (length(stage_line) > 0) sub("!Sample_characteristics_ch1 = age: ", "", stage_line) else NA
  celltype <- if (length(type_line) > 0) sub("!Sample_characteristics_ch1 = cell type: ", "", type_line) else NA
  
  # Return a data frame for each sample (if GSM ID is found)
  if (!is.na(gsm)) {
    return(data.frame(SampleID = gsm, DevelopmentStage = stage, CellType = celltype, stringsAsFactors = FALSE))
  } else {
    return(NULL)
  }
})

# Remove NULL entries and combine all metadata into one data frame
sample_metadata <- sample_metadata[!sapply(sample_metadata, is.null)]
metadata_df <- do.call(rbind, sample_metadata)
# Align metadata with Seurat object columns
metadata_df <- metadata_df[seq_len(ncol(seurat_object_symbols)), ]
# Add extracted metadata as columns to the Seurat object
seurat_object_symbols$sample_id <- metadata_df$SampleID
seurat_object_symbols$development_stage <- metadata_df$DevelopmentStage
seurat_object_symbols$cell_type_soft <- metadata_df$CellType



# Quality Control & Filtering #

# Visualize QC metrics
VlnPlot(seurat_object_symbols, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
FeatureScatter(seurat_object_symbols, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# Filter cells with unusually high RNA counts
seurat_object_symbols <- subset(seurat_object_symbols, subset = nCount_RNA > 400000)
# Print basic summary of the filtered object
summary(seurat_object_symbols)



# Normalization & Feature Selection #

# Identify highly variable genes
seurat_object_symbols <- FindVariableFeatures(seurat_object_symbols, selection.method = "vst")
# Get and plot the top 10 most variable features
top10_features <- head(VariableFeatures(seurat_object_symbols), 10)
plot1 <- VariableFeaturePlot(seurat_object_symbols)
plot2 <- LabelPoints(plot = plot1, points = top10_features, repel = TRUE)



# Scaling Data & PCA #

# Scale the data (center and scale) using the variable genes
seurat_object_symbols <- ScaleData(seurat_object_symbols, features = VariableFeatures(seurat_object_symbols))
# Run PCA on the scaled variable genes
seurat_object_symbols <- RunPCA(seurat_object_symbols, features = VariableFeatures(seurat_object_symbols))
# Visualize PCA loadings and plot PCA reduction
VizDimLoadings(seurat_object_symbols, dims = 1:2, reduction = "pca")
DimPlot(seurat_object_symbols, reduction = "pca") 
# Determine how many PCs to retain
ElbowPlot(seurat_object_symbols)
# Extract PCA embeddings (coordinates for each cell)
pca_embeddings <- Embeddings(seurat_object_symbols, "pca")
# Compute pairwise distances between cells in PCA space
pairwise_distances <- dist(pca_embeddings)



# Unbiased Clustering #

seurat_object_symbols <- FindNeighbors(seurat_object_symbols, dims = 1:20)
seurat_object_symbols <- FindClusters(seurat_object_symbols, resolution = 2.0)
seurat_object_symbols <- RunUMAP(seurat_object_symbols, dims = 1:20)

# Change cluster IDs from 0-9 to 1-10
seurat_object_symbols$seurat_clusters <- as.factor(seurat_object_symbols$seurat_clusters)
levels(seurat_object_symbols$seurat_clusters) <- seq(1, length(levels(seurat_object_symbols$seurat_clusters)))

# Check the new cluster labels
head(seurat_object_symbols$seurat_clusters)



# Identifying Cluster Specific Top 5 Marker Genes in unbiased clusters #

# Find markers for each cell type
seurat_markers <- FindAllMarkers(
  seurat_object_symbols,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

# Extract top 5 markers per cell type
top_markers <- seurat_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5, with_ties = FALSE)

top_marker_genes <- top_markers$gene
FeaturePlot(seurat_object_symbols, features = top_marker_genes) # Not useful with larger number of markers such as 5 for each
seurat_object_symbols <- ScaleData(seurat_object_symbols, features = top_marker_genes)

# Scale and plot heatmap
DoHeatmap(seurat_object_symbols, features = top_marker_genes, group.by = "seurat_clusters") +
  ggtitle("Expression of Top 5 Marker Genes per Unbiased Cluster") +
  theme(axis.text.y = element_text(size = 10))

# View top markers in a tabular format
print(top_marker_genes)




# Identifying Cluster Specific Top 5 Marker Genes in Biased Clusters (Cell types) #

# Set cell type as identity
Idents(seurat_object_symbols) <- "cell_type_soft"  # Replace with your actual cell type column if different

# Find markers for each cell type
celltype_markers <- FindAllMarkers(
  seurat_object_symbols,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

# Extract top 5 markers per cell type
top_celltype_markers <- celltype_markers %>%
  group_by(cluster) %>%  # 'cluster' reflects cell types here
  slice_max(order_by = avg_log2FC, n = 5, with_ties = FALSE)

top_celltype_genes <- top_celltype_markers$gene

# Scale and plot heatmap
seurat_object_symbols <- ScaleData(seurat_object_symbols, features = top_celltype_genes)

DoHeatmap(seurat_object_symbols, features = top_celltype_genes, group.by = "cell_type_soft") +
  ggtitle("Top 5 Marker Genes per Biased Clusters") +
  theme(axis.text.y = element_text(size = 10))

print(top_celltype_genes)




# Adult vs Fetal Comparisons #

# Define three development groups Adult, Fetal_Replicating, Fetal_Quiescent 
seurat_object_symbols$development_stage_simple <- ifelse(
  grepl("fetal.*replicating", seurat_object_symbols$cell_type_soft, ignore.case = TRUE), "Fetal_Replicating",
  ifelse(grepl("fetal.*quiescent", seurat_object_symbols$cell_type_soft, ignore.case = TRUE), "Fetal_Quiescent",
         ifelse(grepl("postnatal", seurat_object_symbols$development_stage, ignore.case = TRUE), "Adult", NA))
)

# Check how many cells fall into each group
table(seurat_object_symbols$development_stage_simple)

# Set new identity class
Idents(seurat_object_symbols) <- "development_stage_simple"

# Run pairwise DE comparisons 

# Fetal_Replicating vs Adult
markers_rep_vs_adult <- FindMarkers(
  seurat_object_symbols,
  ident.1 = "Fetal_Replicating",
  ident.2 = "Adult",
  min.pct = 0.25,
  logfc.threshold = 0.25
)

# Fetal_Quiescent vs Adult
markers_quiescent_vs_adult <- FindMarkers(
  seurat_object_symbols,
  ident.1 = "Fetal_Quiescent",
  ident.2 = "Adult",
  min.pct = 0.25,
  logfc.threshold = 0.25
)

# Fetal_Replicating vs Fetal_Quiescent
markers_rep_vs_quiescent <- FindMarkers(
  seurat_object_symbols,
  ident.1 = "Fetal_Replicating",
  ident.2 = "Fetal_Quiescent",
  min.pct = 0.25,
  logfc.threshold = 0.25
)

# Extract top genes from each comparison (top 20) 

top_genes_rep_vs_adult <- markers_rep_vs_adult %>%
  arrange(p_val) %>%
  filter(!is.na(avg_log2FC)) %>%
  slice_head(n = 20) %>%
  rownames()

top_genes_quiescent_vs_adult <- markers_quiescent_vs_adult %>%
  arrange(p_val) %>%
  filter(!is.na(avg_log2FC)) %>%
  slice_head(n = 20) %>%
  rownames()

top_genes_rep_vs_quiescent <- markers_rep_vs_quiescent %>%
  arrange(p_val) %>%
  filter(!is.na(avg_log2FC)) %>%
  slice_head(n = 20) %>%
  rownames()

# Combine gene list 

top_genes_all <- unique(c(
  top_genes_rep_vs_adult,
  top_genes_quiescent_vs_adult,
  top_genes_rep_vs_quiescent
))

#  print number of genes
cat("Total unique genes selected:", length(top_genes_all), "\n")

# Scale these genes in the Seurat object
seurat_object_symbols <- ScaleData(seurat_object_symbols, features = top_genes_all)
  

# Create a combined group variable (Developmental Stage + Cell Type)
seurat_object_symbols$development_and_celltype <- paste(
  seurat_object_symbols$development_stage_simple,
  seurat_object_symbols$cell_type_soft,
  sep = "_"
)


# Create conditional group labels: only append cell type name to Adult
seurat_object_symbols$development_and_celltype <- ifelse(
  seurat_object_symbols$development_stage_simple == "Adult",
  paste("Adult", seurat_object_symbols$cell_type_soft, sep = "_"),
  seurat_object_symbols$development_stage_simple
)

# Heatmap

DoHeatmap(seurat_object_symbols, features = top_genes_all, group.by = "development_and_celltype") +
  ggtitle("Top Differentially Expressed Genes Between Fetal and Adult Cells") +
  theme(
    axis.text.y = element_text(size = 10), # y-axis gene labels
    plot.title = element_text(size = 15)
  )

print(top_genes_all, group.by = "development_and_celltype")

# Dotplot

DotPlot(
  seurat_object_symbols,
  features = top_genes_all,
  group.by = "development_and_celltype"
) +
  ggtitle("Top Differentially Expressed Genes Between Fetal and Adult Cells") +
  RotatedAxis() +
  theme(
    axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 15)
  )

# Select a smaller set of genes for the boxplot
selected_genes <- top_genes_all[1:5]  # Limit to top 5 genes for clarity


# UMAPs

p_cluster <- DimPlot(seurat_object_symbols, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
  ggtitle("(A) Unbiased Clusters")
p_stage <- DimPlot(seurat_object_symbols, reduction = "umap", group.by = "development_stage_simple", label = TRUE) +
  ggtitle("(B) Developmental Stage")
p_celltype <- DimPlot(seurat_object_symbols, reduction = "umap", group.by = "cell_type_soft", label = TRUE) +
  ggtitle("(C) Biased Clusters")

# Combine horizontally
p_cluster | p_stage | p_celltype +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

# proportion info
prop.table(table(seurat_object_symbols$seurat_clusters, seurat_object_symbols$cell_type_soft), 1)



# Identifying MHC Expression #

DotPlot(seurat_object_symbols, 
        features = c("HLA-A", "HLA-B", "HLA-C", "TAPBP", "CALR", "ERAP1", "B2M", "PDIA3", "HSPA5", "SEC61-A1", "SEC61-A2", "SEC61-B", "SEC61-G"), 
        group.by = "cell_type_soft") +
  ggtitle("MHC1 Gene Expression by Cell Type") +
  RotatedAxis() +
  theme(axis.text.x = element_text(size = 8)) +
  xlab("MHC1 Pathway Genes")  # Custom x-axis label


