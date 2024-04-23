#'
#'clear environment
rm(list = ls())

#' Load the required libraries
library(BiocManager)
library(Biostrings)
library(BiocGenerics)
library(GenomeInfoDb)
library(S4Vectors)
library(stats4)
library(DESeq2)
#'
#' for heatmap
library(gplots)
library(pheatmap)
#'
#' for hierarchical clustering
library(cluster)
#'
library(genefilter)
#'
#' Import dataset into R
rna_data <- read.table(file.choose(), header = TRUE, sep = "\t")

#' 1. Filter for top 1000 most variable genes
gene_expr <- rna_data[, -1]  #' Exclude the 'Gene' column
gene_vars <- apply(gene_expr, 1, var)
top_genes <- rna_data[order(gene_vars, decreasing = TRUE), ][1:min(1000, nrow(rna_data)), ]

#' 2. Compute correlation matrix for top genes
cor_mat <- cor(top_genes[, -1])

#' 3. Perform hierarchical clustering
dist_mat <- as.dist(1 - cor_mat)  # Convert correlation to distance
hclust_res <- hclust(dist_mat, method = "complete")

#' 4. Reorder rows and columns based on clustering
ordered_mat <- cor_mat[hclust_res$order, hclust_res$order]

#' 5. Create heatmap with hierarchical clustering
library(pheatmap)
pheatmap(
  ordered_mat,
  fontsize_row = 8,
  fontsize_col = 8,
  margins = c(5,5),
  main = "Hierarchical Clustering Heatmap (Top 1000 Genes)",
  cluster_rows = TRUE,
  cluster_cols = TRUE
)

