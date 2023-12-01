# load packages ----
library(tidyverse)
library(cowplot)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(Seurat)
library(glmGamPoi)
library(SeuratWrappers)
library(SeuratDisk)
library(reticulate)
options(Seurat.object.assay.version = "v5")
options(future.globals.maxSize = 1e10)

# re-cluster cells from 2580 within seu2
seu2 <- readRDS(file = "RDSfiles/seu2.RDS")
seu2580 <- subset(seu2, subset = orig.ident == "2580")
seu2580 <- NormalizeData(seu2580, verbose = FALSE)
seu2580 <- FindVariableFeatures(seu2580, verbose = FALSE)
seu2580 <- ScaleData(seu2580, verbose = FALSE)
seu2580 <- RunPCA(seu2580, npcs = 10, verbose = FALSE)
seu2580 <- FindNeighbors(seu2580, dims = 1:10, verbose = FALSE)
seu2580 <- FindClusters(seu2580, resolution = 1, verbose = FALSE)
seu2580 <- RunUMAP(seu2580, dims = 1:10, verbose = FALSE)
DimPlot(seu2580, label = TRUE, repel = TRUE) + NoAxes()
DimPlot(seu2580, group.by = "orig.ident") + NoAxes()

FeaturePlot(seu2580,features = "Epcam", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu2580,features = "Ptprc", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu2580,features = "Dcn", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu2580,features = "Muc5ac", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu2580,features = "Muc6", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu2580,features = "Muc4", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu2580,features = "Gif", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu2580,features = "Atp4b", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu2580,features = "Mki67", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()

seu2580 <- RenameIdents(seu2580, 
                     `0` = "Pit0", `1` = "Pit1", `2` = "Pit2", 
                     `3` = "Neck1", `4` = "Ki67", `5` = "Pit3", 
                     `6` = "Pit4", `7` = "Prietal", `8` = "Pit5", `9` = "Neck2", 
                     `10` = "Pit6", `11` = "Chief", `12` ="Pit7")
seu2580$celltype <- Idents(seu2580)
DimPlot(seu2580, label = TRUE, repel = TRUE) + NoAxes()

# save metadata table
seu2580$barcode <- colnames(seu2580)
seu2580$UMAP_1 <- seu2580@reductions$umap@cell.embeddings[,1]
seu2580$UMAP_2 <- seu2580@reductions$umap@cell.embeddings[,2]
write.csv(seu2580@meta.data, file='out/seu2580_metadata.csv', quote=F, row.names=F)

# write expression counts matrix
library(Matrix)
counts_matrix <- LayerData(seu2580, assay='RNA', layer='counts')
writeMM(counts_matrix, file='out/seu2580_counts.mtx')

# write dimesnionality reduction matrix, in this example case pca matrix
write.csv(seu2580@reductions$pca@cell.embeddings, file='out/seu2580_pca.csv', quote=F, row.names=F)

# write gene names
write.table(
  data.frame('gene'=rownames(counts_matrix)),file='out/seu2580_gene_names.csv',
  quote=F,row.names=F,col.names=F
)

saveRDS(seu2580, file = "RDSfiles/seu2580.RDS")
