####################
#### load packages ####
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(Seurat)
library(SeuratWrappers)
options(Seurat.object.assay.version = "v5")


####################
#### load data ####
seu <- readRDS(file = "RDSfiles/seu_01_mt25.RDS")
seu <- JoinLayers(seu)


####################
#### SCTransform clustering without integration ####
seu <- SCTransform(seu, vars.to.regress = c("percent.mt"))
seu <- RunPCA(seu, npcs = 30, verbose = FALSE)
seu <- FindNeighbors(seu, dims = 1:30, verbose = FALSE)
seu <- FindClusters(seu, resolution = 1, verbose = FALSE)
seu <- RunUMAP(seu, dims = 1:30, verbose = FALSE)
DimPlot(seu, label = TRUE, repel = TRUE) + NoAxes()
ggsave("cluster.png", path = "plots/02_sct_wo_integration", width = 5, height = 5, units = "in", dpi = 150)
DimPlot(seu, group.by = "orig.ident") + NoAxes()
ggsave("id.png", path = "plots/02_sct_wo_integration", width = 5, height = 5, units = "in", dpi = 150)


####################
#### feature plots ####
files <- list.files(path = "gene_set/annotation/", full.names = TRUE)
features <- lapply(files, FUN = readLines) %>% unlist()
for(i in 1:length(features)){
  p <- FeaturePlot(seu, features = features[i], cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
  ggsave(paste0(as.character(features[i]), ".png"), plot = p, path = "plots/02_sct_wo_integration", width = 5, height = 5, units = "in", dpi = 150)
}


saveRDS(seu, file = "RDSfiles/seu_03_sct.RDS")