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

# first cluster without integration ----
seu2 <- readRDS(file = "RDSfiles/seu2.RDS")
seu3 <- subset(seu2, subset = gem_prep == "batch1", invert = TRUE)
seu3[["RNA"]]$data <- NULL
seu3[["RNA"]]$scale.data <- NULL
seu3 <- DietSeurat(seu3)
seu3 <- NormalizeData(seu3, verbose = FALSE)
seu3 <- FindVariableFeatures(seu3, verbose = FALSE)
seu3 <- ScaleData(seu3, verbose = FALSE)
seu3 <- RunPCA(seu3, npcs = 10, verbose = FALSE)

seu3 <- FindNeighbors(seu3, dims = 1:10, verbose = FALSE)
seu3 <- FindClusters(seu3, resolution = 0.5, cverbose = FALSE)
seu3 <- RunUMAP(seu3, dims = 1:10, verbose = FALSE)
DimPlot(seu3, label = TRUE, repel = TRUE) + NoAxes()
DimPlot(seu3, group.by = "orig.ident") + NoAxes()
DimPlot(seu3, split.by = "orig.ident", ncol = 2) + NoAxes()
FeaturePlot(seu3,features = "Mki67", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu3,features = "Muc5ac", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu3,features = "Muc6", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu3,features = "Muc4", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu3,features = "Gif", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu3,features = "Atp4b", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()

## split by batch (GEM preparation), sctransformaation with harmony integration ----
## tried dim = 10 or 20
seu2 <- readRDS(file = "RDSfiles/seu2.RDS")
seu3 <- subset(seu2, subset = gem_prep == "batch1", invert = TRUE)
seu3[["RNA"]]$data <- NULL
seu3[["RNA"]]$scale.data <- NULL
seu3 <- DietSeurat(seu3)
seu3[["RNA"]] <- split(seu3[["RNA"]], f = seu3$gem_prep) 

seu3_1 <- SCTransform(seu3, vars.to.regress = "percent.mt")
seu3_1 <- RunPCA(seu3_1, npcs = 10, verbose = FALSE)
seu3_1 <- IntegrateLayers(
  object = seu3_1, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)
seu3_1 <- FindNeighbors(seu3_1, reduction = "harmony", dims = 1:10, verbose = FALSE)
seu3_1 <- FindClusters(seu3_1, resolution = 0.5, cluster.name = "harmony_clusters", verbose = FALSE)
seu3_1 <- RunUMAP(seu3_1, reduction = "harmony", dims = 1:10, verbose = FALSE, reduction.name = "umap.harmony")
DimPlot(seu3_1, reduction = "umap.harmony", label = TRUE, repel = TRUE) + NoAxes()
DimPlot(seu3_1, group.by = "orig.ident") + NoAxes()
DimPlot(seu3_1, split.by = "orig.ident", ncol = 2) + NoAxes()

DefaultAssay(seu3_1) <- "SCT"

FeaturePlot(seu3_1,features = "Epcam", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu3_1,features = "Ptprc", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu3_1,features = "Dcn", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu3_1,features = "Muc5ac", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu3_1,features = "Muc6", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu3_1,features = "Muc4", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu3_1,features = "Gif", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu3_1,features = "Atp4b", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu3_1,features = "Mki67", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()

## split by batch (GEM preparation), lognormalization with harmony integration ----
## tried dim = 10 or 20
seu2 <- readRDS(file = "RDSfiles/seu2.RDS")
seu3 <- subset(seu2, subset = gem_prep == "batch1", invert = TRUE)
seu3[["RNA"]]$data <- NULL
seu3[["RNA"]]$scale.data <- NULL
seu3 <- DietSeurat(seu3)
seu3[["RNA"]] <- split(seu3[["RNA"]], f = seu3$gem_prep)
seu3_2 <- NormalizeData(seu3, verbose = FALSE)
seu3_2 <- FindVariableFeatures(seu3_2, verbose = FALSE)
seu3_2 <- ScaleData(seu3_2, verbose = FALSE)
seu3_2 <- RunPCA(seu3_2, npcs = 10, verbose = FALSE)
seu3_2 <- IntegrateLayers(
  object = seu3_2, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)
seu3_2 <- FindNeighbors(seu3_2, reduction = "harmony", dims = 1:10, verbose = FALSE)
seu3_2 <- FindClusters(seu3_2, resolution = 0.5, cluster.name = "harmony_clusters", verbose = FALSE)
seu3_2 <- RunUMAP(seu3_2, reduction = "harmony", dims = 1:10, verbose = FALSE, reduction.name = "umap.harmony")
DimPlot(seu3_2, reduction = "umap.harmony", label = TRUE, repel = TRUE) + NoAxes()
DimPlot(seu3_2, group.by = "orig.ident") + NoAxes()
DimPlot(seu3_2, split.by = "orig.ident", ncol = 2) + NoAxes()

DefaultAssay(seu3_2) <- "RNA"

FeaturePlot(seu3_2,features = "Epcam", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu3_2,features = "Ptprc", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu3_2,features = "Dcn", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu3_2,features = "Muc5ac", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu3_2,features = "Muc6", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu3_2,features = "Muc4", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu3_2,features = "Gif", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu3_2,features = "Atp4b", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu3_2,features = "Mki67", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()

## lognormalize with scvi integration  ----
# scvi integration did not work with latest versions (Seurat 5.0.1, SeuratObject 5.0.1)
# need to have previous versions as below
# https://github.com/satijalab/seurat/issues/7944
seu2 <- readRDS(file = "RDSfiles/seu2.RDS")
seu3 <- subset(seu2, subset = gem_prep == "batch1", invert = TRUE)
seu3[["RNA"]]$data <- NULL
seu3[["RNA"]]$scale.data <- NULL
seu3 <- DietSeurat(seu3)
seu3[["RNA"]] <- split(seu3[["RNA"]], f = seu3$gem_prep)
seu3_3 <- NormalizeData(seu3, verbose = FALSE)
seu3_3 <- FindVariableFeatures(seu3_3, verbose = FALSE)
seu3_3 <- ScaleData(seu3_3, verbose = FALSE)
seu3_3 <- RunPCA(seu3_3, npcs = 20, verbose = FALSE)
seu3_3 <- IntegrateLayers(
  object = seu3_3, method = scVIIntegration,
  new.reduction = "integrated.scvi",
  conda_env = "/Users/kinoshitahiroto/miniconda3/envs/scvi-env",
  verbose = TRUE
)
seu3_3 <- FindNeighbors(seu3_3, reduction = "integrated.scvi", dims = 1:20, verbose = FALSE)
seu3_3 <- FindClusters(seu3_3, resolution = 0.5, cluster.name = "scvi_clusters", verbose = FALSE)
seu3_3 <- RunUMAP(seu3_3, reduction = "integrated.scvi", dims = 1:20, reduction.name = "umap.scvi", verbose = FALSE)
DimPlot(seu3_3, label = TRUE, repel = TRUE) + NoAxes()
DimPlot(seu3_3, group.by = "orig.ident") + NoAxes()
DimPlot(seu3_3, split.by = "orig.ident", ncol = 2) + NoAxes()

FeaturePlot(seu3_3,features = "Epcam", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu3_3,features = "Ptprc", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu3_3,features = "Dcn", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu3_3,features = "Muc5ac", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu3_3,features = "Muc6", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu3_3,features = "Muc4", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu3_3,features = "Gif", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu3_3,features = "Atp4b", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu3_3,features = "Mki67", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()


#### proceed with lognormalization with harmony integration, dim 10
## rename idents ----
Idents(seu3_2) <- "seurat_clusters"
seu3_2 <- RenameIdents(seu3_2, 
                       `0` = "Pit0", `1` = "Pit1", `2` = "Pit2", 
                       `3` = "Neck1", `4` = "Parietal", `5` = "Ki67", 
                       `6` = "Neck2", `7` = "Pit3", `8` = "Pit4", `9` = "Chief", 
                       `10` = "Immune", `11` = "Parietal", `12` ="Neck2", `13` = "Pit5")
seu3_2$celltype <- Idents(seu3_2)
seu3_2$celltype <- factor(seu3_2$celltype, levels = c("Ki67", "Pit0", "Pit1", "Pit2", "Pit3", "Pit4", "Pit5", "Neck1", "Neck2", "Chief", "Parietal", "Immune"))
Idents(seu3_2) <- "celltype"
DimPlot(seu3_2, label = TRUE, repel = TRUE) + NoAxes()
seu3_2 <- JoinLayers(seu3_2)
saveRDS(seu3_2, file = "RDSfiles/seu3_2.RDS")

# selecting cells by using monocle3 function choose_cells()
seu3_2 <- subset(seu3_2, subset = celltype %in% c("Immune", "Parietal"), invert = TRUE)
seu3_2@reductions$umap <- seu3_2@reductions$umap.harmony
cds <- as.cell_data_set(seu3_2)
cds <- choose_cells(cds)
seu3_2 <- seu3_2[,cds@colData@rownames]
DimPlot(seu3_2, label = TRUE, repel = TRUE) + NoAxes()
seu3_2$orig.ident <- factor(seu3_2$orig.ident, levels = c("1","2","3","4","2547","2580"))
DimPlot(seu3_2, split.by = "orig.ident", ncol = 2) + NoAxes()

FeaturePlot(seu3_2,features = "Mki67", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu3_2,features = "Muc5ac", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu3_2,features = "Muc6", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu3_2,features = "Muc4", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu3_2,features = "Gif", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu3_2,features = "Stmn1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu3_2,features = "Cd44", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu3_2,features = "Gkn3", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()

saveRDS(seu3_2, file = "RDSfiles/seu3_2_selected.RDS")

# monocle3 ----
# https://rpubs.com/mahima_bose/Seurat_and_Monocle3_p
cds <- as.cell_data_set(seu3_2)
fData(cds)$gene_short_name <- rownames(fData(cds))

# assign partitions
recreate.partitions <- c(rep(1, length(cds@colData@rownames)))
names(recreate.partitions) <- cds@colData@rownames
recreate.partitions <- as.factor(recreate.partitions)
cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partitions

# assign clusters
list.cluster <- seu3_2@active.ident
cds@clusters@listData[["UMAP"]][["clusters"]] <- list.cluster

# assign UMAP coordinates
cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- seu3_2@reductions$umap@cell.embeddings

# learn trajectory
cds <- learn_graph(cds, use_partition = F)
plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = F,
           label_branch_points = F, label_roots = T, label_leaves = F,
           group_label_size = 5, graph_label_size = 0, cell_size = 1) + NoAxes()

# order cells in pseudotime
cds <- order_cells(cds, reduction_method = "UMAP", root_cells = colnames(cds[, clusters(cds) == "Ki67"]))
plot_cells(cds, color_cells_by = "pseudotime", label_groups_by_cluster = T,
           label_branch_points = F, label_roots = F, label_leaves = F, cell_size = 1) + NoAxes()

# cds$monocle3_pseudotime <- pseudotime(cds)
# data.pseudo <- as.data.frame(colData(cds))
# ggplot(data.pseudo, aes(monocle3_pseudotime, celltype, fill = celltype)) + geom_boxplot()
# ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(celltype, monocle3_pseudotime), fill = celltype)) + geom_boxplot()
# 
# # find genes that change as a function of pseudotime
# deg <- graph_test(cds, neighbor_graph = "principal_graph")
# deg %>% arrange(q_value) %>% filter(status == "OK") %>% head()
# FeaturePlot(seu3_2, features = c("Cox5b", "Fn1", "Ptma", "Pigr"))

# add pseudotime values to the seurat object
seu3_2$pseudotime <- pseudotime(cds)
FeaturePlot(seu3_2, features = "pseudotime") + NoAxes()
# RidgePlot(seu3_2, features = c("pseudotime"), sort = T, idents = c("Ki67", "Pit10", "Pit3", "Neck1", "Chief"))
# my_genes <- row.names(subset(fData(cds), gene_short_name %in% c("Mki67", "Muc5ac", "Pgc"))) 
# cds_subset <- cds[my_genes,]
# plot_genes_in_pseudotime(cds_subset, color_cells_by = "monocle3_pseudotime" )

# CytoTRACE
counts_matrix <- LayerData(seu3_2, assay='RNA', layer='counts') %>% as.data.frame()
obj_cell_type_anno <- as.data.frame(seu3_2@meta.data$celltype)
results <- CytoTRACE(counts_matrix, ncores = 4)
pheno <- as.character(seu3_2@meta.data$celltype)
names(pheno) <- colnames(seu3_2)
plotCytoTRACE(results, phenotype = pheno)
# Get umap embeddings from Seurat...
# plotCytoTRACE(results, phenotype = pheno, gene = "TIMP1", emb = umapEmb)

## convert Seurat object to anndata manually following the tutorial below ----
# https://smorabit.github.io/tutorials/8_velocyto/

# save metadata table
seu3_2$barcode <- colnames(seu3_2)
seu3_2$UMAP_1 <- seu3_2@reductions$umap@cell.embeddings[,1]
seu3_2$UMAP_2 <- seu3_2@reductions$umap@cell.embeddings[,2]
write.csv(seu3_2@meta.data, file='out/seu3_2_metadata.csv', quote=F, row.names=F)

# write expression counts matrix
library(Matrix)
counts_matrix <- LayerData(seu3_2, assay='RNA', layer='counts')
writeMM(counts_matrix, file='out/seu3_2_counts.mtx')

# write dimesnionality reduction matrix, in this example case pca matrix
write.csv(seu3_2@reductions$pca@cell.embeddings, file='out/seu3_2_pca.csv', quote=F, row.names=F)

# write gene names
write.table(
  data.frame('gene'=rownames(counts_matrix)),file='out/seu3_2_gene_names.csv',
  quote=F,row.names=F,col.names=F
)


