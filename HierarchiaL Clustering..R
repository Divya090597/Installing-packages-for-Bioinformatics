#
# clear environment
rm(list = ls())

# Load the required libraries
library(BiocManager)
library(Biostrings)
library(BiocGenerics)
library(GenomeInfoDb)
library(S4Vectors)
library(stats4)
library(DESeq2)
#
# for heatmap
library(gplots)
library(pheatmap)
#
# for hierarchical clustering
library(cluster)
#
library(genefilter)
#

# Import dataset into R
rna_data <- read.table(file.choose(), header = TRUE, sep = "\t")
rna_data

# Filter for top 1000 most variable genes
gene_expr <- rna_data[, -1]  # Exclude the 'Gene' column
gene_vars <- apply(gene_expr, 1, var)
top_genes <- rna_data[order(gene_vars, decreasing = TRUE), ][1:min(1000, nrow(rna_data)), ]

# Log2 transform the data
log2_data <- log2(top_genes[, -1] + 1)

# Z-score normalize
zdataset <- apply(log2_data, 1, scale)

# Invert Z-score normalization (if needed)
zdataset <- apply(zdataset, 1, rev)

# Convert to matrix and set column names
colnames(zdataset) <- colnames(log2_data)
dataset <- as.matrix(zdataset)

# Compute correlation matrix for top genes
cor_mat <- cor(dataset)

# Perform hierarchical clustering
dist_mat <- as.dist(1 - cor_mat)  # Convert correlation to distance
hclust_res <- hclust(dist_mat, method = "complete")

# Reorder rows and columns based on clustering
ordered_mat <- cor_mat[hclust_res$order, hclust_res$order]

# Create heatmap with hierarchical clustering
library(pheatmap)
pheatmap(
  ordered_mat,
  fontsize_row = 8,
  fontsize_col = 8,
  margins = c(5, 5),
  main = "Hierarchical Clustering Heatmap (Top 1000 Genes)",
  cluster_rows = TRUE,
  cluster_cols = TRUE
)






