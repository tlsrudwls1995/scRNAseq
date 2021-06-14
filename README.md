# scRNAseq_ps47.1
library(dplyr)
library(Seurat)
library(patchwork)

# Load the sample_ps47.1 dataset
sample3.data <- Read10X(data.dir = "../data/sample_ps47.1/sample352/outs/filtered_feature_bc_matrix/")
# Initialize the Seurat object with the raw (non-normalized data).
sample3 <- CreateSeuratObject(counts = sample3.data, project = "pbmc3k", min.cells = 3, min.features = 200)
sample3
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
sample3[["percent.mt"]] <- PercentageFeatureSet(sample3, pattern = "^mt-")
# Visualize QC metrics as a violin plot
VlnPlot(sample3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(sample3, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(sample3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
sample3 <- subset(sample3, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
sample3 <- NormalizeData(sample3)
sample3 <- FindVariableFeatures(sample3, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(sample3), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(sample3)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
all.genes <- rownames(sample3)
sample3 <- ScaleData(sample3, features = all.genes)
sample3 <- RunPCA(sample3, features = VariableFeatures(object = sample3))
# Examine and visualize PCA results a few different ways
print(sample3[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(sample3, dims = 1:2, reduction = "pca")
DimPlot(sample3, reduction = "pca")
DimHeatmap(sample3, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(sample3, dims = 1:15, cells = 500, balanced = TRUE)
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
sample3 <- JackStraw(sample3, num.replicate = 100)
sample3 <- ScoreJackStraw(sample3, dims = 1:20)
JackStrawPlot(sample3, dims = 1:15)
ElbowPlot(sample3)

#Cluster cells
sample3 <- FindNeighbors(sample3, dims = 1:10)
sample3 <- FindClusters(sample3, resolution = 0.5)
# Look at cluster IDs of the first 5 cells
head(Idents(sample3), 5)
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
sample3 <- RunUMAP(sample3, dims = 1:10)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(sample3, reduction = "umap")
saveRDS(sample3, file = "../data/sample_ps47.1/sample3.rds")

# Find differentially expressed features(cluster biomarkers)
# find all markers of cluster 2
cluster2.markers <- FindMarkers(sample3, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(sample3, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)
# find markers for every cluster compared to all remaining cells, report only the positive ones
pbmc.markers <- FindAllMarkers(sample3, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
cluster0.markers <- FindMarkers(sample3, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
VlnPlot(sample3, features = c("Arl6ip1", "Ube2c"))
# you can plot raw counts as well
VlnPlot(sample3, features = c("Crip2", "Lgals1"), slot = "counts", log = TRUE)
FeaturePlot(sample3, features = c("Actb", "Dbi", "Tnnc1", "Myl7", "Actc1", "Rbp1", "Hba-a1", "Hbb-y",
                                  "Apela"))
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(sample3, features = top10$gene) + NoLegend()
new.cluster.ids <- c("MS cells", "'Cluster1'", "VMSC", "'Cluster3'", "'Cluster4'", "'Cluster5'")
names(new.cluster.ids) <- levels(sample3)
sample3 <- RenameIdents(sample3, new.cluster.ids)
DimPlot(sample3, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
saveRDS(sample3, file = "../data/sample_ps47.1/sample3_final.rds")

