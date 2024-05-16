d

#' Import dataset into R
rna_data <- read.table(file.choose(), header = TRUE, sep = "\t")

#' Additional Preprocessing
# Assuming 'exp' is your dataset
dataset <- rna_data

# Log2 transformation
dataset <- log2(dataset + 1)

# Z-score normalization
zdataset <- apply(dataset, 1, scale)

# Reverse order within each row
zdataset <- apply(zdataset, 1, rev)

# Set column names
colnames(zdataset) <- colnames(dataset)

#' 1. Filter for top 1000 most variable genes
gene_expr <- zdataset[, -1]  #' Exclude the 'Gene' column
gene_vars <- apply(gene_expr, 1, var)
top_genes <- zdataset[order(gene_vars, decreasing = TRUE), ][1:min(1000, nrow(zdataset)), ]

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
