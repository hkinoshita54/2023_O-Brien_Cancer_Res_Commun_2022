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
#### SCTransform clustering with harmony ####
seu[["RNA"]] <- split(seu[["RNA"]], f = seu$gem_prep)
seu <- NormalizeData(seu, verbose = FALSE)
seu <- FindVariableFeatures(seu, verbose = FALSE)
seu <- ScaleData(seu, verbose = FALSE)
seu <- RunPCA(seu, npcs = 20, verbose = FALSE)
seu <- IntegrateLayers(
  object = seu, method = HarmonyIntegration,
  orig.reduction = "pca", 
  new.reduction = "harmony")
seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:20, verbose = FALSE)
seu <- FindClusters(seu, reduction = "harmony", resolution = 1, verbose = FALSE)
seu <- RunUMAP(seu, reduction = "harmony", dims = 1:20, verbose = FALSE)
DimPlot(seu, label = TRUE, repel = TRUE) + NoAxes()
ggsave("cluster.png", path = "plots/04_lognorm_harmony", width = 5, height = 5, units = "in", dpi = 150)
DimPlot(seu, group.by = "orig.ident") + NoAxes()
ggsave("id.png", path = "plots/04_lognorm_harmony", width = 5, height = 5, units = "in", dpi = 150)


####################
#### feature plots ####
files <- list.files(path = "gene_set/annotation/", full.names = TRUE)
features <- lapply(files, FUN = readLines) %>% unlist()
for(i in 1:length(features)){
  p <- FeaturePlot(seu, features = features[i], cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
  ggsave(paste0(as.character(features[i]), ".png"), plot = p, path = "plots/04_lognorm_harmony", width = 5, height = 5, units = "in", dpi = 150)
}


saveRDS(seu, file = "RDSfiles/seu_05_lognorm_harmony.RDS")