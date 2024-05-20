
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

# Filter for top 1000 most variable genes using MAD
# Calculate the MAD for each gene, log2 transformed. Select the top 1000 genes based on MAD.

gene_expr <- rna_data[, -1]  # Exclude the 'Gene' column
mad_values <- apply(log2(gene_expr + 1), 1, mad)
top_genes_index <- order(mad_values, decreasing = TRUE)[1:min(1000, length(mad_values))]
top_genes <- rna_data[top_genes_index, ]

# Log2 transform the data
log2_data <- log2(top_genes[, -1] + 1)

# Z-score normalize and invert Z-score normalization (if needed)
zdataset <- t(apply(log2_data, 1, scale))
zdataset <- t(apply(zdataset, 1, rev))

# Convert to matrix and set column names
colnames(zdataset) <- colnames(log2_data)
dataset <- as.matrix(zdataset)
dataset[is.na(dataset)] <- 0  # Replace NA values with 0
dataset <- dataset[rowSums(dataset != 0) > 0, ]  # Remove rows with all zeros

# Perform hierarchical clustering
dist_mat <- dist(dataset)  # Use Euclidean distance
hclust_res <- hclust(dist_mat, method = "complete")

# Reorder rows and columns based on clustering
ordered_mat <- dataset[hclust_res$order, ]

# Create heatmap with hierarchical clustering
pheatmap(
  ordered_mat,
  fontsize_row = 5,
  fontsize_col = 5,
  margins = c(5, 5),
  main = "Hierarchical Clustering Heatmap (Top 1000 Genes)",
  cluster_rows = TRUE,
  cluster_cols = TRUE
)






