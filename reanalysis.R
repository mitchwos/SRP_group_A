library(lattice)
library(flexmix)
library(scde)
library(GEOquery)
library(SingleCellExperiment)

# read gene matrix csv
df <- read.csv('/home/mjw85/Documents/SRP/Group_Project/gene_count_matrix.csv', 
               header = TRUE, row.names=1)
# remove extra info off column names
colnames(df) <- sub('_.*', '', colnames(df))
# order columns 
order = sort(colnames(df))
df <- df[,order]

# get metadata files from GEO
gse <- getGEO('GSE67835', GSEMatrix=TRUE)
# access files in gse
gse1 <- gse[[1]]
gse2 <- gse[[2]]
# extract metadata from files into dataframes
meta1 <- pData(phenoData(gse1))
meta2 <- pData(phenoData(gse2))
# combine dataframes
meta <- rbind(meta1, meta2)

# extract and rename relevant columns
meta$age <- meta$`age:ch1`
meta$cell_type <- meta$`cell type:ch1`
# discard all other columns 
meta <- meta[c('age','cell_type')]
# order rows
new_order = sort(rownames(meta))
meta <- meta[new_order,]
# remove years and weeks from age values
meta$age <- sub('postnatal \\d+ years?', 'postnatal', meta$age)
meta$age <- sub('prenatal 16-18 W', 'prenatal', meta$age)
# make age and cell_type factors
meta$age <- factor(meta$age, levels=c('postnatal','prenatal'))
meta$cell_type <- factor(meta$cell_type, levels=c('oligodendrocytes','hybrid',
                                                    'OPC','astrocytes','neurons',
                                                    'microglia','endothelial',
                                                    'fetal_replicating',
                                                    'fetal_quiescent'))
# construct sce object
sce <- SingleCellExperiment(assays=list(counts=df), colData=meta)

# list of genes with postnatal value
post_cells <- row.names(meta[meta[['age']]=='postnatal',])
# list of postnatal and prenatal matching df columns
col_lab <- ifelse(names(df) %in% post_cells, 'postnatal','prenatal')
# convert df to matrix
mat <- data.matrix(df)
# add attribute to matrix of age factor 
attr(mat,'age') <- factor(col_lab, levels=c('postnatal','prenatal'))
age <- attr(mat,'age')
# save as file
save(mat, file='/home/mjw85/Documents/SRP/Group_Project/mat.RData')
# build error model with age groups
em <- scde.error.models((counts=mat), groups=age, n.cores=1,
                        threshold.segmentation=TRUE, save.crossfit.plots=FALSE,
                        save.model.plots=FALSE, verbose=1)
# save as file
save(em, file='/home/mjw85/Documents/SRP/Group_Project/em.RData')

# parallel and scde
library(parallel)
require(scde)
# weighted distance measure method
m.prior <- scde.expression.prior(models=em, counts=mat)
jp <- scde.posteriors(models=em, mat, m.prior, 
                      return.individual.posterior.modes=TRUE, n.cores=1)
jp$jp.modes <- log(as.numeric(colnames(jp$jp)))[max.col(jp$jp)]
p.mode.fail <- scde.failure.probability(models=em, magnitudes=jp$jp.modes)
p.self.fail <- scde.failure.probability(models=em, counts=mat)
# weight matrix
matw <- 1-sqrt(p.self.fail*sqrt(p.self.fail*p.mode.fail))
# magnitude matrix (using individual posterior modes)
mat <- log10(exp(jp$modes)+1);
# weighted distance
cell.names <- colnames(mat); names(cell.names) <- cell.names;
require(boot)
mode.fail.dist <- as.dist(1-do.call(rbind,mclapply(cell.names,function(nam1) {
  unlist(lapply(cell.names,function(nam2) {
    corr(cbind(mat[,nam1],mat[,nam2]),w=sqrt(sqrt(matw[,nam1]*matw[,nam2])))
  }))
},mc.cores=1)),upper=F);

# show dendrogram
clust <- hclust(mode.fail.dist, method='ward.D2');
plot(as.dendrogram(clust))

# Rtsne
library(Rtsne)
tsne <- Rtsne(mode.fail.dist, dims=3, perplexity=30)

# BIC
library(mclust)
BIC <- mclustBIC(tsne$Y, G=1:40)
# save as file
save(BIC, file='/home/mjw85/Documents/SRP/Group_Project/BIC.RData')
# Gaussian mixture model 
gmm <- Mclust(tsne$Y, G=1:40)
# save as file
save(gmm, file='/home/mjw85/Documents/SRP/Group_Project/gmm.RData')
# plot BIC and add optimal cluster number
plot(BIC)
abline(v = gmm$G)

# scatterplot3d
library(scatterplot3d)
# clusters from gmm 
clusters <- gmm$classification
# cluster df 
c_df <- data.frame(tsne$Y)
c_df$cluster <- clusters
# make clusters factors
c_df$cluster <- factor(c_df$cluster, levels=c('1','2','3','4','5','6','7','8','9','10'))
# assign clusters colors
colors <- c('#f8453d','#f8913d','#c6c42b','#74c30d','#1eb84a','#33e5d7','#1c97eb',
            '#ac4be1','#e742ea','#f24dad')
colors <- colors[as.numeric(c_df$cluster)]
# 3D scatter plot of clusters
sc3d <- scatterplot3d(c_df[-4], pch=20, color=colors, box=FALSE, xlab='', ylab='',
              zlab='')
legend(sc3d$xyz.convert(11, 17, 7), legend=levels(c_df$cluster), 
       col=c('#f8453d','#f8913d','#c6c42b','#74c30d','#1eb84a','#33e5d7','#1c97eb',
             '#ac4be1','#e742ea','#f24dad'), pch=20)

# boxplot of uncertainty by cluster
df_list <- list(cluster=gmm$classification, uncertainty=gmm$uncertainty)
uc_df <- as.data.frame(df_list)

x1 <- uc_df$uncertainty[uc_df$cluster==1]
x2 <- uc_df$uncertainty[uc_df$cluster==2]
x3 <- uc_df$uncertainty[uc_df$cluster==3]
x4 <- uc_df$uncertainty[uc_df$cluster==4]
x5 <- uc_df$uncertainty[uc_df$cluster==5]
x6 <- uc_df$uncertainty[uc_df$cluster==6]
x7 <- uc_df$uncertainty[uc_df$cluster==7]
x8 <- uc_df$uncertainty[uc_df$cluster==8]
x9 <- uc_df$uncertainty[uc_df$cluster==9]
x10 <- uc_df$uncertainty[uc_df$cluster==10]

boxplot(x1, x2, x3, x4, x5, x6, x7, x8, x9, x10)

# check number of cells in each cluster
cluster_numbers <- data.frame(table(c_df$cluster))
colnames(cluster_numbers) <- c('cluster','number of cells')
# save as csv
write.csv(cluster_numbers,'/home/mjw85/Documents/SRP/Group_Project/cluster_nos.csv',
          row.names=FALSE)

# FactoMineR
library(FactoMineR)
library(factoextra)
# PCA cluster plot
res.pca <- PCA(c_df[-4], graph=FALSE)
fviz_pca_ind(res.pca, label='none', habillage=c_df$cluster)

# add tsne to sce 
mtsne <- data.matrix(tsne$Y)
reducedDim(sce, 'tsne') <- mtsne
# add clusters to metadata
colData(sce)$clusters <- as.factor(clusters)

# empty list
cluster_deg <- list()
# function to apply to each cluster
diff_exp <- function(x) {
  # factor groups of cluster and all other clusters
  cg <- factor(ifelse(colData(sce)$clusters==x, x, 'other'))
  names(cg) <- colnames(mat)
  # prior expression model
  prior <- scde.expression.prior(models=em, counts=mat)
  # expression difference
  ediff <- scde.expression.difference(em, mat, prior, groups=cg,
                                      n.randomizations=100, n.cores=1, 
                                      verbose=1)
  # convert z-score to p-values
  p.value <- 2*pnorm(abs(ediff$Z), lower.tail=F)
  # convert corrected z-score to adjusted p-values
  p.adjust <- 2*pnorm(abs(ediff$cZ), lower.tail=F)
  # find and order most significant genes
  sig.genes <- which(p.adjust<0.05)
  ord <- order(p.adjust[sig.genes])
  # create ordered matrix
  de <- cbind(ediff[sig.genes,1:3],p.adjust[sig.genes])[ord,]
  # assign column names
  colnames(de) <- c('Lower bound','log2 fold change','Upper bound','p-value')
  # store results globally
  cluster_deg[[x]] <<- de
}
# apply function to all clusters
invisible(lapply(1:10, diff_exp))
# save as file
save(cluster_deg, file='/home/mjw85/Documents/SRP/Group_Project/cluster_DEG.RData')

# add DEG by clusters to sce object
metadata(sce)$cluster_DEG <- cluster_deg
# save as file
save(sce, file='/home/mjw85/Documents/SRP/Group_Project/sce.RData')

# empty list
genes_list <- list()
# function to extract top 20 DEGs from each cluster
top_20 <- function(x){
  genes <- row.names(cluster_deg[[x]])[1:20]
  genes_list[[x]] <<- genes
}
# apply function to all clusters
lapply(1:10, top_20)
# df of clusters' top 20 DEGs
clusters_20 <- data.frame(genes_list)
colnames(clusters_20) <- c(1:10)
# save as csv
write.csv(clusters_20,'/home/mjw85/Documents/SRP/Group_Project/top_20.csv',
          row.names=FALSE)

# empty list
cell_deg <- list()
# function to apply to each cell type
cell_diff_exp <- function(x) {
  # factor groups of cell type and all other cell types
  cells <- factor(ifelse(colData(sce)$cell_type==x, x, 'other'))
  names(cells) <- colnames(mat)
  # prior expression model
  prior <- scde.expression.prior(models=em, counts=mat)
  # expression difference
  ediff <- scde.expression.difference(em, mat, prior, groups=cells,
                                      n.randomizations=100, n.cores=1, 
                                      verbose=1)
  # convert z-score to p-values
  p.value <- 2*pnorm(abs(ediff$Z), lower.tail=F)
  # convert corrected z-score to adjusted p-values
  p.adjust <- 2*pnorm(abs(ediff$cZ), lower.tail=F)
  # find and order most significant genes
  sig.genes <- which(p.adjust<0.05)
  ord <- order(p.adjust[sig.genes])
  # create ordered matrix
  de <- cbind(ediff[sig.genes,1:3],p.adjust[sig.genes])[ord,]
  # assign column names
  colnames(de) <- c('Lower bound','log2 fold change','Upper bound','p-value')
  # store results globally
  cell_deg[[x]] <<- de
}
# list of cell types
cell_types <- c('oligodendrocytes','hybrid','OPC','astrocytes','neurons',
'microglia','endothelial','fetal_replicating','fetal_quiescent')
# apply function to all clusters
invisible(lapply(cell_types, cell_diff_exp))
# save as file
save(cell_deg, file='/home/mjw85/Documents/SRP/Group_Project/cells_DEG.RData')

# add DEG by cell type to sce object
metadata(sce)$cell_DEG <- cell_deg
# save as file
save(sce, file='/home/mjw85/Documents/SRP/Group_Project/sce.RData')

# empty list
genes_list2 <- list()
# function to extract top 100 DEGs from each cell type
top_100 <- function(x){
  genes <- row.names(cell_deg[[x]])[1:100]
  genes_list2[[x]] <<- genes
}
# apply function to all cell types
lapply(cell_types, top_100)
# df of cell types' top 100 DEGs
cells_100 <- data.frame(genes_list2)
colnames(cells_100) <- cell_types
# save as csv
write.csv(cells_100,'/home/mjw85/Documents/SRP/Group_Project/cells_100.csv',
          row.names=FALSE)

# packages
library(tidyr)
library(dplyr)
# lengthen data with cell type and gene columns
cell_df <- pivot_longer(cells_100, cols=everything(), names_to='cell type',
                        values_to='gene')
# lengthen data with cluster and gene columns
cluster_df <- pivot_longer(clusters_20, cols=everything(), names_to='cluster',
                           values_to='gene')
# match genes in cells_100 to genes in cluster_20
matches <- cell_df %>%
  filter(gene %in% cluster_df$gene) %>%
  inner_join(cluster_df, by='gene', relationship='many-to-many')
# save as file
save(matches, file='/home/mjw85/Documents/SRP/Group_Project/matches.RData')

# roperators
library(roperators)
# read gene marker csv
markers <- read.tsv('/home/mjw85/Documents/SRP/Group_Project/gene_markers.csv')
# function to clean gene names so they will match
clean_gene <- function(x) {
  x %>%
    as.character() %>%
    trimws() %>%
    toupper() %>%
    gsub("[^A-Z0-9]", "", .)
}
# clean cluster gene names
cluster_df <- cluster_df %>%
  mutate(gene=clean_gene(gene))
# clean marker gene names
markers <- markers %>%
  mutate(gene=clean_gene(gene))
# match genes in cluster_df and markers
marker_matches <- markers %>%
  filter(gene %in% cluster_df$gene) %>%
  inner_join(cluster_df, by='gene', relationship='many-to-many')
# save as file
save(marker_matches, file='/home/mjw85/Documents/SRP/Group_Project/marker_matches.RData' )

# count matrix
count_mat <- counts(sce)
# scuttle
library(scuttle)
# apply log normalization to sce
sce <- scuttle::logNormCounts(sce)
# adds assay of logcounts to sce object
assay(sce, 'logcounts') <- logcounts(sce)

# markers genes from paper's heat map
mouse_genes <- c('DAAM2','ASPA','MAL','SEC14L5','MAP6D1','DPYD','PPP1R14A',
                 'GJB1','FA2H','MAG','CDK18','LGI3','SHC4','UGT8','KLK6',
                 'KCNH8','GPR37','MOBP','LPAR1','ADAMTS4','ERMN','OPALIN',
                 'CLDN11','PLEKHB1','GSN','GRM3','CNP','MBP','PLP1','SLC14A1',
                 'GLIS3','GLI3','PPP1R3C','CHRDL1','CYBRD1','CTH','SORCS2',
                 'ITGB4','RNF43','NWD1','PAQR6','C16orf89','ALDH1L1','TRIM66',
                 'HGF','CBS','ITGA7','SLC30A10','SLC4A4','FGFR3','BMPR1B',
                 'ATP13A4','AQP4','GPR183','CCL4','CD83','LAPTM5','CSF1R',
                 'HLA-DRA','BCL2A1','CD14','CCL2','APOLD1','TM4SF1','FLT1',
                 'A2M','PDGFRA','LHFPL3','MEGF11','PCDH15','KCNK1','KIAA1324',
                 'LNX1','NELL1','COBL','SLITRK1','DPYSL5','C14orf37','DLX1',
                 'DLX6','GLRA2','DLX2','DLX5','SLC10A4','EGFR','SST','PNOC',
                 'NXPH1','BCL11A','DCN','TMEM130','CNTN4','CDO1','NFASC',
                 'LRRTM3','GRIA3','RELN')
# clean gene names
mouse_genes <- toupper(trimws(sub('\\..*','', mouse_genes)))
# remove fetal cells from sce object
sce_adult <- sce[, !colData(sce)$cell_type %in% c('fetal_replicating',
                                                  'fetal_quiescent')]
# log-normalized expression data of adult cells
mouse_exp <- assay(sce_adult,'logcounts')
# clean gene names
row.names(mouse_exp) <- toupper(trimws(sub('\\..*', '', row.names(mouse_exp))))
# only keep marker genes
mouse_exp <- mouse_exp[row.names(mouse_exp) %in% mouse_genes, ]
# scale matrix
scaled_mouse_exp <- t(scale(t(mouse_exp)))
# save as file
save(scaled_mouse_exp, file='/home/mjw85/Documents/SRP/Group_Project/select_genes.RData')
# distance matrix
cell_dist <- dist(t(scaled_mouse_exp))
# hierarchical clustering
cell_hclust <- hclust(cell_dist, method='ward.D2')
# plot dendrogram
plot(cell_hclust, labels = FALSE)

# tibble
library(tibble)
# tibble of cells and 7 biased clusters
biased_clusters <- cutree(cell_hclust, k=7) %>%
  enframe() %>%
  dplyr::rename(cell=name, cluster=value)
# add biased clusters to sce adult
colData(sce_adult)$biased_clusters <- factor(biased_clusters$cluster)
# save as file
save(sce_adult, file='/home/mjw85/Documents/SRP/Group_Project/sce_adult.RData')

# annotation of biased and unbiased clusters for heat map
cluster_anno <- data.frame(biased = (colData(sce_adult)$biased_clusters),
                           unbiased = (colData(sce_adult)$clusters))
# naming row names after cells 
row.names(cluster_anno) <- row.names(colData(sce_adult))
# save as csv
gene_anno <- read.csv('/home/mjw85/Documents/SRP/Group_Project/cell_marker_genes.csv', 
                       header = TRUE, row.names=1)
# order genes in expression matrix to match gene annotation AKA by cell type
scaled_mouse_exp <- scaled_mouse_exp[rownames(gene_anno), ]
# list for factor levels
cluster_order <- c('1','2','3','4','5','6','7')
# order cells by clusters
ordered_cells <- rownames(cluster_anno)[
  order(factor(cluster_anno$unbiased, levels=cluster_order))
]
# apply order to expression matrix columns
scaled_mouse_exp <- scaled_mouse_exp[, ordered_cells]
# apply order to cluster annotation rows
cluster_anno <- cluster_anno[ordered_cells, , drop=FALSE]
# assign colors to marker gene cell types
gene_colors <- list(
  cell.type = c('oligodendrocyte genes'='#fa9e5d','astrocyte genes'='#f3ee6e',
                'microglia genes'='#65ff54','endothelial genes'='#f7b1df',
                'OPC genes'='#8bafff','neuron genes'='#ff7e6a')
)
# convert biased clusters to characters and factorize
cluster_anno$biased <- as.factor(as.character(cluster_anno$biased))
# convert unbiased clusters to characters and factorize
cluster_anno$unbiased <- as.factor(as.character(cluster_anno$unbiased))
# assign colors to clusters
cluster_colors <- list(
  biased = c('1'='#e35e49','2'='#f09e48','3'='#f1ea66','4'='#8fe357',
             '5'='#3a79d3','6'='#be6fea','7'='#ffa8dc'),
  unbiased = c('1'='#e35e49','2'='#f09e48','3'='#3a79d3','4'='#f1ea66',
               '5'='#8fe357','6'='#be6fea','7'='#ffa8dc')
)
# pheatmap
library(pheatmap)
# create heat map 
pheatmap(scaled_mouse_exp,
         annotation_row=gene_anno,
         annotation_colors=c(gene_colors,cluster_colors),
         annotation_col=cluster_anno,
         cluster_rows=FALSE,
         cluster_cols=FALSE,
         show_rownames=TRUE,
         show_colnames=FALSE,
         annotation_legend=FALSE,
         annotation_names_row=FALSE)
