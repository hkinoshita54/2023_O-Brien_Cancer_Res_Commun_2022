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

## split by batch (GEM preparation), sctransformaation with harmony integration ----
## tried dim = 10 or 20
seu2 <- readRDS(file = "RDSfiles/seu2.RDS")
seu2[["RNA"]]$data <- NULL
seu2[["RNA"]]$scale.data <- NULL
seu2 <- DietSeurat(seu2)
seu2[["RNA"]] <- split(seu2[["RNA"]], f = seu2$gem_prep) 
seu2_1 <- SCTransform(seu2, vars.to.regress = "percent.mt")
seu2_1 <- RunPCA(seu2_1, npcs = 20, verbose = FALSE)
seu2_1 <- IntegrateLayers(
  object = seu2_1, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)
seu2_1 <- FindNeighbors(seu2_1, reduction = "harmony", dims = 1:20, verbose = FALSE)
seu2_1 <- FindClusters(seu2_1, resolution = 0.5, cluster.name = "harmony_clusters", verbose = FALSE)
seu2_1 <- RunUMAP(seu2_1, reduction = "harmony", dims = 1:20, verbose = FALSE, reduction.name = "umap.harmony")
DimPlot(seu2_1, reduction = "umap.harmony", label = TRUE, repel = TRUE) + NoAxes()
DimPlot(seu2_1, group.by = "orig.ident") + NoAxes()
DimPlot(seu2_1, split.by = "orig.ident", ncol = 2) + NoAxes()

DefaultAssay(seu2_1) <- "SCT"

FeaturePlot(seu2_1,features = "Epcam", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu2_1,features = "Ptprc", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu2_1,features = "Dcn", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu2_1,features = "Muc5ac", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu2_1,features = "Muc6", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu2_1,features = "Muc4", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu2_1,features = "Gif", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu2_1,features = "Atp4b", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu2_1,features = "Mki67", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()

## split by batch (GEM preparation), lognormalization with harmony integration ----
## tried dim = 10 or 20
seu2 <- readRDS(file = "RDSfiles/seu2.RDS")
seu2[["RNA"]]$data <- NULL
seu2[["RNA"]]$scale.data <- NULL
seu2 <- DietSeurat(seu2)
seu2[["RNA"]] <- split(seu2[["RNA"]], f = seu2$gem_prep)
seu2_2 <- NormalizeData(seu2, verbose = FALSE)
seu2_2 <- FindVariableFeatures(seu2_2, verbose = FALSE)
seu2_2 <- ScaleData(seu2_2, verbose = FALSE)
seu2_2 <- RunPCA(seu2_2, npcs = 20, verbose = FALSE)
seu2_2 <- IntegrateLayers(
  object = seu2_2, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)
seu2_2 <- FindNeighbors(seu2_2, reduction = "harmony", dims = 1:20, verbose = FALSE)
seu2_2 <- FindClusters(seu2_2, resolution = 2, cluster.name = "harmony_clusters", verbose = FALSE)
seu2_2 <- RunUMAP(seu2_2, reduction = "harmony", dims = 1:20, verbose = FALSE, reduction.name = "umap.harmony")
DimPlot(seu2_2, reduction = "umap.harmony", label = TRUE, repel = TRUE) + NoAxes()
DimPlot(seu2_2, group.by = "orig.ident") + NoAxes()
DimPlot(seu2_2, split.by = "orig.ident", ncol = 2) + NoAxes()

DefaultAssay(seu2_2) <- "RNA"

FeaturePlot(seu2_2,features = "Epcam", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu2_2,features = "Ptprc", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu2_2,features = "Dcn", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu2_2,features = "Muc5ac", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu2_2,features = "Muc6", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu2_2,features = "Muc4", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu2_2,features = "Gif", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu2_2,features = "Atp4b", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu2_2,features = "Mki67", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()

## split by batch (GEM preparation), lognormalize with scvi integration  ----
# scvi integration did not work with latest versions (Seurat 5.0.1, SeuratObject 5.0.1)
# need to have previous versions as below
# https://github.com/satijalab/seurat/issues/7944
## tried dim = 20　and 30
seu2 <- readRDS(file = "RDSfiles/seu2.RDS")
seu2[["RNA"]]$data <- NULL
seu2[["RNA"]]$scale.data <- NULL
seu2 <- DietSeurat(seu2)
seu2[["RNA"]] <- split(seu2[["RNA"]], f = seu2$gem_prep)
seu2_3 <- NormalizeData(seu2, verbose = FALSE)
seu2_3 <- FindVariableFeatures(seu2_3, verbose = FALSE)
seu2_3 <- ScaleData(seu2_3, verbose = FALSE)
seu2_3 <- RunPCA(seu2_3, npcs = 30, verbose = FALSE)
seu2_3 <- IntegrateLayers(
  object = seu2_3, method = scVIIntegration,
  new.reduction = "integrated.scvi",
  conda_env = "/Users/kinoshitahiroto/miniconda3/envs/scvi-env",
  verbose = TRUE
)
seu2_3 <- FindNeighbors(seu2_3, reduction = "integrated.scvi", dims = 1:30, verbose = FALSE)
seu2_3 <- FindClusters(seu2_3, resolution = 0.5, cluster.name = "scvi_clusters", verbose = FALSE)
seu2_3 <- RunUMAP(seu2_3, reduction = "integrated.scvi", dims = 1:30, reduction.name = "umap.scvi", verbose = FALSE)
DimPlot(seu2_3, label = TRUE, repel = TRUE) + NoAxes()
DimPlot(seu2_3, group.by = "orig.ident") + NoAxes()
DimPlot(seu2_3, split.by = "orig.ident", ncol = 2) + NoAxes()

FeaturePlot(seu2_3,features = "Epcam", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu2_3,features = "Ptprc", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu2_3,features = "Dcn", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu2_3,features = "Muc5ac", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu2_3,features = "Muc6", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu2_3,features = "Muc4", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu2_3,features = "Gif", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu2_3,features = "Atp4b", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu2_3,features = "Mki67", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()


#### proceed with lognorm_harmony_dim20 ----
## rename idents ----
seu2_2 <- RenameIdents(seu2_2, 
                     `0` = "Pit4", `1` = "Pit2", `2` = "Pit4", `3` = "Neck3", `4` = "Pit5", 
                     `5` = "Ki67", `6` = "Pit5", `7` = "Pit1", `8` = "Parietal", `9` = "Neck2", 
                     `10` = "Neck1", `11` = "Ki67", `12` ="Uncomitted", `13` = "Pit5", `14` = "Chief", 
                     `15` = "Ki67", `16` = "Pit4", `17` = "Pit6", `18` = "Pit4", `19` = "Parietal", `20` = "Parietal",
                     `21` = "?", `22` = "Pit5", `23` = "Pit3", `24` = "Parietal", `25` = "Pit2")
seu2_2$celltype <- Idents(seu2_2)
seu2_2$celltype <- factor(seu2_2$celltype, levels = c("Ki67", "Pit1", "Pit2", "Pit3", "Pit4", "Pit5", "Pit6", "Neck1", "Neck2", "Neck3", "Chief", "Parietal", "Uncomitted", "?"))
Idents(seu2_2) <- seu2_2$celltype
DimPlot(seu2_2, label = TRUE, repel = TRUE) + NoAxes()
DimPlot(seu2_2, split.by = "orig.ident", ncol = 2) + NoAxes()
saveRDS(seu2_2, file = "RDSfiles/seu2_2.RDS")

# selecting cells by using monocle3 function choose_cells() ----
seu2_2 <- readRDS(file = "RDSfiles/seu2_2.RDS")
seu2_2 <- JoinLayers(seu2_2)
seu2_2 <- subset(seu2_2, subset = celltype %in% c("Parietal", "Uncomitted", "?"), invert = TRUE)
seu2_2@reductions$umap <- seu2_2@reductions$umap.harmony
library(monocle3)
cds <- as.cell_data_set(seu2_2)
cds <- choose_cells(cds)
seu2_2 <- seu2_2[,cds@colData@rownames]
DimPlot(seu2_2, label = TRUE, repel = TRUE) + NoAxes()
seu2_2$orig.ident <- factor(seu2_2$orig.ident, levels = c("2266","2204","1","2","3","4","2547","2580"))
DimPlot(seu2_2, split.by = "orig.ident", ncol = 2) + NoAxes()

FeaturePlot(seu2_2,features = "Muc5ac", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu2_2,features = "Muc6", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu2_2,features = "Muc4", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu2_2,features = "Gif", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu2_2,features = "Mki67", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu2_2,features = "Stmn1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu2_2,features = "Cd44", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()

saveRDS(seu2_2, file = "RDSfiles/seu2_2_selected.RDS")

## monocle3 ----
# https://rpubs.com/mahima_bose/Seurat_and_Monocle3_p
cds <- as.cell_data_set(seu2_2)
fData(cds)$gene_short_name <- rownames(fData(cds))

# assign partitions
recreate.partitions <- c(rep(1, length(cds@colData@rownames)))
names(recreate.partitions) <- cds@colData@rownames
recreate.partitions <- as.factor(recreate.partitions)
cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partitions

# assign clusters
list.cluster <- seu2_2@active.ident
cds@clusters@listData[["UMAP"]][["clusters"]] <- list.cluster

# assign UMAP coordinates
cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- seu2_2@reductions$umap@cell.embeddings

# learn trajectory
cds <- learn_graph(cds, use_partition = F)
plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = F,
           label_branch_points = F, label_roots = T, label_leaves = F,
           group_label_size = 5, graph_label_size = 0, cell_size = 1) + NoAxes()

# order cells in pseudotime
cds <- order_cells(cds, reduction_method = "UMAP", root_cells = colnames(cds[, clusters(cds) == "Ki67"]))
plot_cells(cds, color_cells_by = "pseudotime", label_groups_by_cluster = T,
           label_branch_points = F, label_roots = F, label_leaves = F, cell_size = 1) + NoAxes()



####

## convert Seurat object to anndata manually following the tutorial below ----
# https://smorabit.github.io/tutorials/8_velocyto/

# save metadata table
seu2_2$barcode <- colnames(seu2_2)
seu2_2$UMAP_1 <- seu2_2@reductions$umap@cell.embeddings[,1]
seu2_2$UMAP_2 <- seu2_2@reductions$umap@cell.embeddings[,2]
write.csv(seu2_2@meta.data, file='out/seu2_2_metadata.csv', quote=F, row.names=F)

# write expression counts matrix
library(Matrix)
counts_matrix <- LayerData(seu2_2, assay='RNA', layer='counts')
writeMM(counts_matrix, file='out/seu2_2_counts.mtx')

# write dimesnionality reduction matrix, in this example case pca matrix
write.csv(seu2_2@reductions$pca@cell.embeddings, file='out/seu2_2_pca.csv', quote=F, row.names=F)

# write gene names
write.table(
  data.frame('gene'=rownames(counts_matrix)),file='out/seu2_2_gene_names.csv',
  quote=F,row.names=F,col.names=F
)


