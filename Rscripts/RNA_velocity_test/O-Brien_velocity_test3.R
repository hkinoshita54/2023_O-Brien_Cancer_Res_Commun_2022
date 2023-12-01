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
seu_l <- as.Seurat(x = ldat)
seu_l$orig.ident <- "1"

# change cell names according to cellranger output
colnames(seu_l) <- sapply(strsplit(colnames(seu_l), ":"), "[", 2)
colnames(seu_l) <- gsub("x", "-1", colnames(seu_l))

# load counts from cellranger output
mtx <- dir(path = "./data/GSE224840_RAW/", pattern = paste0("matrix", 1, "\\.mtx"), full.names = TRUE)
features <- dir(path = "./data/GSE224840_RAW/", pattern = paste0("features", 1, "\\.tsv"), full.names = TRUE)
barcodes <- dir(path = "./data/GSE224840_RAW/", pattern = paste0("barcodes", 1, "\\.tsv"), full.names = TRUE)
cts <- ReadMtx(mtx = mtx, features = features, cells = barcodes)

# check all the cell names are included in both, order cell names in cts
all(colnames(cts) %in% colnames(seu_l))
all(colnames(cts) == colnames(seu_l))
cts <- cts[, colnames(seu_l)]
all(colnames(cts) == colnames(seu_l))

# create seurat object from cts, then process as normal
seu <- CreateSeuratObject(counts = cts, project = "1", min.cells = 3, min.features = 200)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^mt-")
# VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
seu <- subset(seu, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 25)
seu <- SCTransform(seu)
seu <- RunPCA(seu)
seu <- RunUMAP(seu, dims = 1:20)
seu <- FindNeighbors(seu, dims = 1:20)
seu <- FindClusters(seu)
DimPlot(seu, label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "Atp4b", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()

# save as anndata, then merge with data from loom file in python (or with anndata from loom file generated below) 
DefaultAssay(seu) <- "RNA"
SaveH5Seurat(seu, filename = "out/seu_test3.h5Seurat", overwrite = TRUE)
Convert("out/seu_test3.h5Seurat", dest = "h5ad", overwrite = TRUE)

# generate anndata from loom file, only selecting cells in seu from cts
seu_l <- seu_l[, colnames(seu)]
seu_l[["RNA"]] <- seu_l[["spliced"]]
DefaultAssay(seu_l) <- "RNA"
SaveH5Seurat(seu_l, filename = "out/seu_l.h5Seurat", overwrite = TRUE)
Convert("out/seu_l.h5Seurat", dest = "h5ad", overwrite = TRUE)
