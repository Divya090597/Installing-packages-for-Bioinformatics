#install.packages
install.packages("gplots")
BiocManager::install("genefilter")


# Load required libraries
library(gplots)  # for heatmap
library(cluster) # for hierarchical clustering
library(genefilter)

# Import dataset into R
rna_data <- read.table(file.choose(), header = TRUE, sep = "\t")
rna_data


# 1. Filter for top 1000 most variable genes
gene_expr <- rna_data[, -1]  # Exclude the 'Gene' column
gene_vars <- apply(gene_expr, 1, var)
top_genes <- rna_data[order(gene_vars, decreasing = TRUE), ][1:1000, ]

# 2. Perform hierarchical clustering
dist_matrix <- dist(top_genes[, -1])  # Exclude the 'Gene' column
hclust_result <- hclust(dist_matrix)

# Convert top_genes to a numeric matrix
top_genes_matrix <- as.matrix(top_genes[,-1])


# Define custom color palette
custom_palette <- colorRampPalette(c("blue", "white", "red"))  # Blue to white to red gradient

# Generate heatmap
heatmap(top_genes_matrix, col = heatmap_col, scale = "row", 
        Rowv = as.dendrogram(hclust_result), Colv = NA, 
        margins = c(5,5), main = "Hierarchical Clustering Heatmap",
        xlab = "Samples", ylab = "Genes",
        cexRow = 0.3, cexCol = 0.3)
          

                                   





