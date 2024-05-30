
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

# Reading the files manually
barcodes_cd36ko <- read.delim("~/Datasets/cd36ko/GSM5220548_cd36ko_barcodes.tsv.gz", header = FALSE)
features_cd36ko <- read.delim("~/Datasets/cd36ko/GSM5220548_cd36ko_features.tsv.gz", header = FALSE)
matrix_cd36ko <- Matrix::readMM("~/Datasets/cd36ko/GSM5220548_cd36ko_matrix.mtx.gz")

barcodes_cd36wt <- read.delim("~/Datasets/cd36wt/GSM5220547_cd36wt_barcodes.tsv.gz", header = FALSE)
features_cd36wt <- read.delim("~/Datasets/cd36wt/GSM5220547_cd36wt_features.tsv.gz", header = FALSE)
matrix_cd36wt <- Matrix::readMM("~/Datasets/cd36wt/GSM5220547_cd36wt_matrix.mtx.gz")

#Step1: Load the required packages
library(Seurat)
library(Matrix)

# Rename the files to match Read10X expectations
file.rename("~/Datasets/cd36wt/GSM5220547_cd36wt_barcodes.tsv.gz", "~/Datasets/cd36wt/barcodes.tsv.gz")
file.rename("~/Datasets/cd36wt/GSM5220547_cd36wt_features.tsv.gz", "~/Datasets/cd36wt/features.tsv.gz")
file.rename("~/Datasets/cd36wt/GSM5220547_cd36wt_matrix.mtx.gz", "~/Datasets/cd36wt/matrix.mtx.gz")

# Import and preprocess the data
# Ensure that the data is in appropriate format such as matrix or data frame, with cells as rows and genes as columns.
# Read the 10X Genomics data
cd36wt <- Read10X(data.dir = "~/Datasets/cd36wt")

# Rename the files to match Read10X expectations
file.rename("~/Datasets/cd36ko/GSM5220548_cd36ko_barcodes.tsv.gz", "~/Datasets/cd36ko/barcodes.tsv.gz")
file.rename("~/Datasets/cd36ko/GSM5220548_cd36ko_features.tsv.gz", "~/Datasets/cd36ko/features.tsv.gz")
file.rename("~/Datasets/cd36ko/GSM5220548_cd36ko_matrix.mtx.gz", "~/Datasets/cd36ko/matrix.mtx.gz")

cd36ko <- Read10X(data.dir = "~/Datasets/cd36ko")

# Create Seurat Objects to store the data and metadata
Cd36_ko <- CreateSeuratObject(counts = cd36ko, project = "cd36ko", min.cells = 3, min.features = 200)

Cd36_wt <- CreateSeuratObject(counts = cd36wt, project = "cd36wt", min.cells = 3, min.features = 200)

#View Seurat objects
Cd36_ko
colnames(Cd36_ko[])
rownames(Cd36_ko[])
view(Cd36_ko)

Cd36_wt
colnames(Cd36_wt[])
rownames(Cd36_wt[])
view(Cd36_wt)
# STANDARD PREPROCESSING WORKFLOW

# Filter out low quality cells and genes and normalize the data
# Calculate mitochondrial QC metrics
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Cd36_ko[["percent.mt"]] <- PercentageFeatureSet(Cd36_ko, pattern = "^MT-")

Cd36_wt[["percent.mt"]] <- PercentageFeatureSet(Cd36_wt, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Cd36_ko, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

VlnPlot(Cd36_wt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1_ko <- FeatureScatter(Cd36_ko, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2_ko <- FeatureScatter(Cd36_ko, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1_ko + plot2_ko

plot1_wt <- FeatureScatter(Cd36_wt, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2_wt <- FeatureScatter(Cd36_wt, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1_wt + plot2_wt

#NORMALIZING THE DATA

Cd36_ko <- NormalizeData(Cd36_ko, normalization.method = "LogNormalize", scale.factor = 10000)

Cd36_wt <- NormalizeData(Cd36_wt, normalization.method = "LogNormalize", scale.factor = 10000)


#IDENTIFICATION OF HIGHLY VARIABLE FEATURES (feature selection)

Cd36_ko <- FindVariableFeatures(Cd36_ko, selection.method = "vst", nfeatures = 2000)

Cd36_wt <- FindVariableFeatures(Cd36_wt, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10_ko <- head(VariableFeatures(Cd36_ko), 10)

top10_wt <- head(VariableFeatures(Cd36_wt), 10)

# plot variable features with and without labels
plot3_ko <- VariableFeaturePlot(Cd36_ko)
plot4_ko <- LabelPoints(plot = plot3_ko, points = top10_ko, repel = TRUE)
plot3_ko + plot4_ko

plot3_wt <- VariableFeaturePlot(Cd36_wt)
plot4_wt <- LabelPoints(plot = plot3_wt, points = top10_wt, repel = TRUE)
plot3_wt + plot4_wt

# Print individual plots
print(plot3_ko)
print(plot4_ko)

print(plot3_wt)
print(plot4_wt)

# Save the combined plot to a PNG file
png(filename = "combined_plot_ko.png", width = 800, height = 600)
print(plot3_ko + plot4_ko)
dev.off()

png(filename = "combined_plot_wt.png", width = 800, height = 600)
print(plot3_wt + plot4_wt)
dev.off()

# SCALING THE DATA

all.genes <- rownames(Cd36_ko)
Cd36_ko <- ScaleData(Cd36_ko, features = all.genes)

all.genes <- rownames(Cd36_wt)
Cd36_wt <- ScaleData(Cd36_wt, features = all.genes)

# PERFORM LINEAR DIMENTIONAL REDUCTION

Cd36_ko <- RunPCA(Cd36_ko, features = VariableFeatures(object = Cd36_ko))

Cd36_wt <- RunPCA(Cd36_wt, features = VariableFeatures(object = Cd36_wt))

# Examine and visualize PCA results a few different ways
#Cd36_ko
print(Cd36_ko[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(Cd36_ko, dims = 1:2, reduction = "pca")

DimPlot(Cd36_ko, reduction = "pca") + NoLegend()

DimHeatmap(Cd36_ko, dims = 1, cells = 500, balanced = TRUE)

DimHeatmap(Cd36_ko, dims = 1:15, cells = 500, balanced = TRUE)

#Cd36_wt
print(Cd36_wt[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(Cd36_wt, dims = 1:2, reduction = "pca")

DimPlot(Cd36_wt, reduction = "pca") + NoLegend()

DimHeatmap(Cd36_wt, dims = 1, cells = 500, balanced = TRUE)

DimHeatmap(Cd36_wt, dims = 1:15, cells = 500, balanced = TRUE)

# DETERMINE THE DIMENTIONALITY OF THE DATASET.

ElbowPlot(Cd36_ko)

ElbowPlot(Cd36_wt)

# CLUSTER THE CELLS

Cd36_ko <- FindNeighbors(Cd36_ko, dims = 1:10)
Cd36_ko <- FindClusters(Cd36_ko, resolution = 0.5)

Cd36_wt <- FindNeighbors(Cd36_wt, dims = 1:10)
Cd36_wt <- FindClusters(Cd36_wt, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(Cd36_ko), 5)

head(Idents(Cd36_wt), 5)

#RUN NON-LINEAR DIMENTIONAL REDUCTION

Cd36_ko <- RunUMAP(Cd36_ko, dims = 1:10)

Cd36_wt <- RunUMAP(Cd36_wt, dims = 1:10)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Cd36_ko, reduction = "umap")

DimPlot(Cd36_wt, reduction = "umap")

#Finding differentially expressed features (cluster biomarkers)
# find all markers of cluster 2
cluster2_ko.markers <- FindMarkers(Cd36_ko, ident.1 = 2)
head(cluster2_ko.markers, n = 5)

cluster2_wt.markers <- FindMarkers(Cd36_wt, ident.1 = 2)
head(cluster2_wt.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5_ko.markers <- FindMarkers(Cd36_ko, ident.1 = 5, ident.2 = c(0, 3))
head(cluster5_ko.markers, n = 5)

cluster5_wt.markers <- FindMarkers(Cd36_wt, ident.1 = 5, ident.2 = c(0, 3))
head(cluster5_wt.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
Cd36_ko.markers <- FindAllMarkers(Cd36_ko, only.pos = TRUE)
Cd36_ko.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)
cluster0_ko.markers <- FindMarkers(Cd36_ko, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

Cd36_wt.markers <- FindAllMarkers(Cd36_wt, only.pos = TRUE)
Cd36_wt.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)
cluster0_wt.markers <- FindMarkers(Cd36_wt, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

VlnPlot(Cd36_ko, features = c("Gzmc", "Ccl24"))
VlnPlot(Cd36_wt, features = c("Gzmc", "Gzmd"))

# you can plot raw counts as well
VlnPlot(Cd36_ko, features = c("Gzmc", "Ccl24"), slot = "counts", log = TRUE)

VlnPlot(Cd36_wt, features = c("Gzmf", "Gzmd"), slot = "counts", log = TRUE)

# List all available features
available_features_ko <- rownames(Cd36_ko)
print(head(available_features_ko))

available_features_wt <- rownames(Cd36_wt)
print(head(available_features_wt))

# List the variable features if any are calculated
variable_features_ko <- VariableFeatures(Cd36_ko)
print(head(variable_features_ko))

variable_features_wt <- VariableFeatures(Cd36_wt)
print(head(variable_features_wt))

FeaturePlot(Cd36_ko, features = c("Rpl12", "Top2a", "Ccl7", "Plin2", "Cd8a", "Cdca3", "Cd2", "Hp", "mt-Co1", "Neat1", "Kmo",
                               "Ccl22"))
FeaturePlot(Cd36_wt, features = c("Bcl2", "Cd68", "Cdca8", "Nkg7", "Hp", "Nop2", "Cd2", "Bcl11a", "Tyr",
                                  "Ccl22"))
Cd36_ko.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10_ko
DoHeatmap(Cd36_ko, features = top10$gene) + NoLegend()

Cd36_wt.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10_wt
DoHeatmap(Cd36_wt, features = top10$gene) + NoLegend()


# Merging the Seurat objects
mergeddata <- merge(Cd36_ko, y = Cd36_wt, add.cell.ids = c("Cd36_ko", "Cd36_wt"), project = "mergeddata")
ls()

#View Seurat objects
mergeddata
colnames(mergeddata[])
rownames(mergeddata[])
View(mergeddata)

View(mergeddata@meta.data)

# STANDARD PREPROCESSING WORKFLOW

# Filter out low quality cells and genes and normalize the data
# Calculate mitochondrial QC metrics

#### QC metrics: "nfeature_RNA", "nCount_RNA", "percent.mt"

# Low quality cells or empty droplets often have very few genes
# Cell doublets or multiplets have high values of nfeature_RNA & nCount_RNA
# Low quality cells often have high percentage of mitochondrial genes

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
# Check the assays present in the object
names(mergeddata@assays)

# Check the structure of the RNA assay
str(mergeddata@assays$RNA)

View(mergeddata@meta.data)

range(mergeddata$nFeature_RNA)

range(mergeddata$nCount_RNA)

mergeddata <- PercentageFeatureSet(mergeddata, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(mergeddata, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot5_ID <- FeatureScatter(mergeddata, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot6_ID <- FeatureScatter(mergeddata, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot5_ID + plot6_ID

mergeddata <- subset(mergeddata, subset = nFeature_RNA > 500 & nFeature_RNA < 8000 & nCount_RNA < 10000)

VlnPlot(mergeddata, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

mergeddata

# Data Normalization

mergeddata <- NormalizeData(mergeddata,
                                normalization.method = "LogNormalize", scale.factor = 10000)
# Identification of highly variable features
mergeddata <- FindVariableFeatures(mergeddata, selection.method = "vst", nfeatures = 2000)

# Top features and featureplot
VariableFeaturePlot(mergeddata)

# Data Scaling
mergeddata <- ScaleData(mergeddata)

all.genes <- rownames(mergeddata)
mergeddata <- ScaleData(mergeddata, features = all.genes)

# Perform PCA
mergeddata <- RunPCA(mergeddata)

# Examine and visualize PCA results a few different ways

DimPlot(mergeddata, reduction = "pca", dims = c(1,2))
DimPlot(mergeddata, reduction = "pca", dims = c(1,10))
DimPlot(mergeddata, reduction = "pca", dims = c(1,50))

DimHeatmap(mergeddata, dims = 1, cells = 500, balanced = TRUE)

DimHeatmap(mergeddata, dims = 1:15, cells = 500, balanced = TRUE)


# Determine the dimensionality of the dataset
mergeddata <- JackStraw(mergeddata, num.replicate = 100)
mergeddata <- ScoreJackStraw(mergeddata, dims = 1:20)
JackStrawPlot(mergeddata, dims = 1:20)

ElbowPlot(mergeddata)


# Cluster the cells
mergeddata <- FindNeighbors(mergeddata, dims = 1:10)
mergeddata <- FindClusters(mergeddata, resolution = 0.1)
mergeddata <- FindClusters(mergeddata, resolution = 0.3)
mergeddata <- FindClusters(mergeddata, resolution = 0.5)

# Run UMAP/tSNE for visualization
mergeddata <- RunUMAP(mergeddata, dims = 1:10)
DimPlot(mergeddata, reduction = "umap", label = TRUE, repel = TRUE)

mergeddata <- RunTSNE(object = mergeddata)
DimPlot(object = mergeddata, reduction = "tsne")


##### Setup Seurat objects
#Split the object into a list of 2 repeats
#mergeddata_list <- SplitObject(mergeddata, split.by = 'orig.ident')

#normalize and identify variable features for each dataset independently
mergeddata_list <- lapply(X = mergeddata_list, FUN = function(x){
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  })

# Select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = list(Cd36_ko,Cd36_wt))

# Find integration anchors
integration_anchors <- FindIntegrationAnchors(object.list = list(Cd36_ko,Cd36_wt), anchor.features = features)

# Perform integration to create an 'integrated' data assay
integrated_data <- IntegrateData(anchorset = integration_anchors, dims = 1:30)

# Perform a integration analysis
DefaultAssay(integrated_data) <- "integrated"

integrated_data <- ScaleData(integrated_data, verbose = FALSE)
integrated_data <- RunPCA(integrated_data, npcs = 30, verbose = FALSE)
integrated_data <- FindNeighbors(integrated_data, reduction = "pca", dims = 1:30)
integrated_data <- FindClusters(integrated_data, resolution = 0.5)
integrated_data <- RunUMAP(integrated_data, reduction = "pca", dims = 1:30)

# Visualization
DimPlot(integrated_data, reduction = "umap", label = TRUE)

# Compare
plot7 <- DimPlot(mergeddata, reduction = "umap", group.by = 'orig.ident')
plot8 <- DimPlot(integrated_data, reduction = "umap", group.by = 'orig.ident')
plot7+plot8


# Find markers for each cluster
integrated.markers <- FindAllMarkers(integrated_data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)


# Check the layers within each assay
print(integrated_data@assays$integrated$data)

FeaturePlot(integrated_data, features = c("Birc5", "C1qb", "Lyz2", "Chaf1a", "Inhba", "Ccl5", "Sardh", "Tnfrsf4", "Prima1",
                                  "Ncr1", "Syn3", "Unc5c"))

# Top 10 markers for each cluster
top10_ID <- integrated.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

# Heatmap of top markers
DoHeatmap(integrated_data, features = top10_ID$gene) + NoLegend()

Idents(integrated_data)

new.cluster.ids <- c("CD4+", "CD8+", "T cell-Mono", "CD8+", "TGF-beta", "NK",
                     "CD4", "CD4+", "CD8+", "NK", "10", "11")
names(new.cluster.ids) <- levels(integrated_data)
integrated_data <- RenameIdents(integrated_data, new.cluster.ids)

# UMAP plot with labels
DimPlot(integrated_data, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

View(integrated_data@meta.data)


#PSEUDO BULKFLOW ANALYSIS

# Load required packages
BiocManager::install("ExperimentHub")

library(ExperimentHub)
library(DESeq2)

# Since the required data is already a seurat object

# Check available assays in the Seurat object
View(integrated_data@assays)

integrated_data_Filtered <- subset(integrated_data, subset = nFeature_RNA > 500 & nFeature_RNA < 8000 & nCount_RNA > 10000 & percent.mt < 5)

VlnPlot(mergeddata, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)


BiocManager::install('multtest')
install.packages('metap')


