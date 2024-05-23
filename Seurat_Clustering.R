

library(dplyr)
library(patchwork)
library(tidyverse)

# Install Seurat if not already installed
install.packages("Seurat")

# Install missing dependencies
install.packages("spatstat.random")
install.packages("spatstat.data")
install.packages("spatstat.geom")
install.packages("spatstat.core")

# Load the Seurat package
library(Seurat)

# Reading the files manually
barcodes_cd36ko <- read.delim("~/Datasets/cd36ko/GSM5220548_cd36ko_barcodes.tsv.gz", header = FALSE)
features_cd36ko <- read.delim("~/Datasets/cd36ko/GSM5220548_cd36ko_features.tsv.gz", header = FALSE)
matrix_cd36ko <- Matrix::readMM("~/Datasets/cd36ko/GSM5220548_cd36ko_matrix.mtx.gz")


# Load the required packages
library(Seurat)
library(Matrix)

# Rename the files to match Read10X expectations
file.rename("~/Datasets/cd36wt/GSM5220547_cd36wt_barcodes.tsv.gz", "~/Datasets/cd36wt/barcodes.tsv.gz")
file.rename("~/Datasets/cd36wt/GSM5220547_cd36wt_features.tsv.gz", "~/Datasets/cd36wt/features.tsv.gz")
file.rename("~/Datasets/cd36wt/GSM5220547_cd36wt_matrix.mtx.gz", "~/Datasets/cd36wt/matrix.mtx.gz")

# Read the 10X Genomics data
cd36wt <- Read10X(data.dir = "~/Datasets/cd36wt")

# Rename the files to match Read10X expectations
file.rename("~/Datasets/cd36ko/GSM5220548_cd36ko_barcodes.tsv.gz", "~/Datasets/cd36ko/barcodes.tsv.gz")
file.rename("~/Datasets/cd36ko/GSM5220548_cd36ko_features.tsv.gz", "~/Datasets/cd36ko/features.tsv.gz")
file.rename("~/Datasets/cd36ko/GSM5220548_cd36ko_matrix.mtx.gz", "~/Datasets/cd36ko/matrix.mtx.gz")

cd36ko <- Read10X(data.dir = "~/Datasets/cd36ko")

# Create Seurat Objects
Cd36_ko <- CreateSeuratObject(counts = cd36ko, project = "cd36ko", min.cells = 3, min.features = 200)

Cd36_wt <- CreateSeuratObject(counts = cd36wt, project = "cd36wt", min.cells = 3, min.features = 200)

#View Seurat objects
Cd36_ko
colnames(Cd36_ko[])
rownames(Cd36_ko[])
view(Cd36_ko)

#STANDARD PREPROCESSING WORKFLOW

# Calculate mitochondrial QC metrics
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Cd36_ko[["percent.mt"]] <- PercentageFeatureSet(Cd36_ko, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Cd36_ko, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(Cd36_ko, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Cd36_ko, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#NORMALIZING THE DATA

Cd36_ko <- NormalizeData(Cd36_ko, normalization.method = "LogNormalize", scale.factor = 10000)
Cd36_ko <- NormalizeData(Cd36_ko)

#IDENTIFICATION OF HIGHLY VARIABLE FEATURES (feature selection)

Cd36_ko <- FindVariableFeatures(Cd36_ko, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(Cd36_ko), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(Cd36_ko)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# Print individual plots
print(plot1)
print(plot2)

# Save the combined plot to a PNG file
png(filename = "combined_plot.png", width = 800, height = 600)
print(plot1 + plot2)
dev.off()

# SCALING THE DATA

all.genes <- rownames(Cd36_ko)
Cd36_ko <- ScaleData(Cd36_ko, features = all.genes)

# PERFORM LINEAR DIMENTIONAL REDUCTION

Cd36_ko <- RunPCA(Cd36_ko, features = VariableFeatures(object = Cd36_ko))

# Examine and visualize PCA results a few different ways
print(Cd36_ko[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(Cd36_ko, dims = 1:2, reduction = "pca")

DimPlot(Cd36_ko, reduction = "pca") + NoLegend()

DimHeatmap(Cd36_ko, dims = 1, cells = 500, balanced = TRUE)

DimHeatmap(Cd36_ko, dims = 1:15, cells = 500, balanced = TRUE)

# DETERMINE THE DIMENTIONALITY OF THE DATASET.

ElbowPlot(Cd36_ko)

# CLUSTER THE CELLS

Cd36_ko <- FindNeighbors(Cd36_ko, dims = 1:10)
Cd36_ko <- FindClusters(Cd36_ko, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(Cd36_ko), 5)

#RUN NON-LINEAR DIMENTIONAL REDUCTION

Cd36_ko <- RunUMAP(Cd36_ko, dims = 1:10)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Cd36_ko, reduction = "umap")

#Finding differentially expressed features (cluster biomarkers)
# find all markers of cluster 2
cluster2.markers <- FindMarkers(Cd36_ko, ident.1 = 2)
head(cluster2.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(Cd36_ko, ident.1 = 5, ident.2 = c(0, 3))
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
Cd36_ko.markers <- FindAllMarkers(Cd36_ko, only.pos = TRUE)
Cd36_ko.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)
cluster0.markers <- FindMarkers(Cd36_ko, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

VlnPlot(Cd36_ko, features = c("Cdca8", "Stmn1"))
# you can plot raw counts as well
VlnPlot(Cd36_ko, features = c("Gzmf", "Ccl1"), slot = "counts", log = TRUE)

# List all available features
available_features <- rownames(Cd36_ko)
print(head(available_features))

# List the variable features if any are calculated
variable_features <- VariableFeatures(Cd36_ko)
print(head(variable_features))

FeaturePlot(Cd36_ko, features = c("Bcl2", "Hcst", "Rps24", "Rps27", "Lck", "Rpl17", "Rpl13a", "Fau",
                               "Gm42418"))
Cd36_ko.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(Cd36_ko, features = top10$gene) + NoLegend()

# ASIGNING CELL TYPE IDENTITY TO CLUSTERS

new.cluster.ids <- c("Ribosomal Protein", "DNA topoisomerase", "Macrophage", "Sebocyte", "NK", "Stem",
                     "T", "Hepatic", "Cytochrome", "lncRNA", "Cyclin", "Neuro")
names(new.cluster.ids) <- levels(Cd36_ko)
Cd36_ko <- RenameIdents(Cd36_ko, new.cluster.ids)
DimPlot(Cd36_ko, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
