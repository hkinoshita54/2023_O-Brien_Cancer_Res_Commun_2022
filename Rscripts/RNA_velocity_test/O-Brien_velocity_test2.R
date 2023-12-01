# load packages ----
library(tidyverse)
library(cowplot)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(Seurat)
library(SeuratWrappers)
library(SeuratDisk)
options(Seurat.object.assay.version = "v3") # needed to run as.Seurat()

# load loom files ----
# https://htmlpreview.github.io/?https://github.com/satijalab/seurat-wrappers/blob/master/docs/scvelo.html
ldat <- ReadVelocity(file = "data/loom/1.loom")
seu_velo <- as.Seurat(x = ldat)
seu_velo[["RNA"]] <- seu_velo[["spliced"]]

# load counts from cellranger output
mtx <- dir(path = "./data/GSE224840_RAW/", pattern = paste0("matrix", 1, "\\.mtx"), full.names = TRUE)
features <- dir(path = "./data/GSE224840_RAW/", pattern = paste0("features", 1, "\\.tsv"), full.names = TRUE)
barcodes <- dir(path = "./data/GSE224840_RAW/", pattern = paste0("barcodes", 1, "\\.tsv"), full.names = TRUE)
cts <- ReadMtx(mtx = mtx, features = features, cells = barcodes)

# change cell names according to cellranger output
colnames(seu_velo) <- sapply(strsplit(colnames(seu_velo), ":"), "[", 2)
colnames(seu_velo) <- gsub("x", "-1", colnames(seu_velo))

# check all the cell names are included in both, order cell names in cts
all(colnames(cts) %in% colnames(seu_velo))
all(colnames(cts) == colnames(seu_velo))
cts <- cts[, colnames(seu_velo)]
all(colnames(cts) == colnames(seu_velo))

# keep features which are in both cts and loom files
seu_velo <- DietSeurat(seu_velo, features = rownames(cts))
cts <- cts[rownames(seu_velo),] 

# create seurat object from cts, then put it to "RNA" of the object
seu <- CreateSeuratObject(counts = cts, project = "1")
seu_velo[["RNA"]] <- seu[["RNA"]]


# clustering as in the vignette
seu_velo <- SCTransform(seu_velo)
seu_velo <- RunPCA(seu_velo)
seu_velo <- RunUMAP(seu_velo, dims = 1:20)
seu_velo <- FindNeighbors(seu_velo, dims = 1:20)
seu_velo <- FindClusters(seu_velo)
DefaultAssay(seu_velo) <- "RNA"
SaveH5Seurat(seu_velo, filename = "out/seu_velo.h5Seurat")
Convert("out/seu_velo.h5Seurat", dest = "h5ad")

DimPlot(seu_velo, label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_velo,features = "Muc5ac", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
