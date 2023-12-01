# load packages ----
library(tidyverse)
library(cowplot)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
options(Seurat.object.assay.version = "v3") # needed to run as.Seurat()

# load loom files ----
# https://htmlpreview.github.io/?https://github.com/satijalab/seurat-wrappers/blob/master/docs/scvelo.html
ldat <- ReadVelocity(file = "data/loom/1.loom")
seu_sp <- as.Seurat(x = ldat) # use "spliced" as "RNA"
seu_sp[["RNA"]] <- seu_sp[["spliced"]]
seu_sp <- SCTransform(seu_sp)
seu_sp <- RunPCA(seu_sp)
seu_sp <- RunUMAP(seu_sp, dims = 1:20)
seu_sp <- FindNeighbors(seu_sp, dims = 1:20)
seu_sp <- FindClusters(seu_sp)
DefaultAssay(seu_sp) <- "RNA"
SaveH5Seurat(seu_sp, filename = "out/seu.h5Seurat")
Convert("out/seu.h5Seurat", dest = "h5ad")

DimPlot(seu_sp, label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_sp,features = "Dcn", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
