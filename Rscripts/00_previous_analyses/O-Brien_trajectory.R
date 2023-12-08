# load packages ----
library(tidyverse)
library(cowplot)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(Seurat)
library(SeuratWrappers)
library(SeuratDisk)
options(Seurat.object.assay.version = "v5")
library(reticulate)
use_condaenv("cytotrace")
library(CytoTRACE)
library(monocle3)
options(future.globals.maxSize = 1e10)

# load data
seu2 <- readRDS(file = "RDSfiles/seu2.RDS")

# selecting cells by using monocle3 function choose_cells()
seu2 <- subset(seu2, subset = celltype %in% c("Uncomitted", "Parietal1", "Parietal2", "Parietal3", "??"), invert = TRUE)
cds <- as.cell_data_set(seu2)
cds <- choose_cells(cds)
seu2 <- seu2[,cds@colData@rownames]
DimPlot(seu2, label = TRUE, repel = TRUE) + NoAxes()
DimPlot(seu2, split.by = "orig.ident", ncol = 2) + NoAxes()

FeaturePlot(seu2,features = "Mki67", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu2,features = "Muc5ac", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu2,features = "Muc6", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu2,features = "Muc4", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu2,features = "Gif", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu2,features = "Stmn1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu2,features = "Cd44", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu2,features = "Gkn3", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()

# monocle3 ----
# https://rpubs.com/mahima_bose/Seurat_and_Monocle3_p
cds <- as.cell_data_set(seu2)
fData(cds)$gene_short_name <- rownames(fData(cds))

# assign partitions
recreate.partitions <- c(rep(1, length(cds@colData@rownames)))
names(recreate.partitions) <- cds@colData@rownames
recreate.partitions <- as.factor(recreate.partitions)
cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partitions

# assign clusters
list.cluster <- seu2@active.ident
cds@clusters@listData[["UMAP"]][["clusters"]] <- list.cluster

# assign UMAP coordinates
cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- seu2@reductions$umap@cell.embeddings

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
# FeaturePlot(seu2, features = c("Cox5b", "Fn1", "Ptma", "Pigr"))

# add pseudotime values to the seurat object
seu2$pseudotime <- pseudotime(cds)
FeaturePlot(seu2, features = "pseudotime") + NoAxes()
# RidgePlot(seu2, features = c("pseudotime"), sort = T, idents = c("Ki67", "Pit10", "Pit3", "Neck1", "Chief"))
# my_genes <- row.names(subset(fData(cds), gene_short_name %in% c("Mki67", "Muc5ac", "Pgc"))) 
# cds_subset <- cds[my_genes,]
# plot_genes_in_pseudotime(cds_subset, color_cells_by = "monocle3_pseudotime" )

# CytoTRACE
counts_matrix <- LayerData(seu2, assay='RNA', layer='counts') %>% as.data.frame()
obj_cell_type_anno <- as.data.frame(seu2@meta.data$celltype)
results <- CytoTRACE(counts_matrix, ncores = 4)
pheno <- as.character(seu2@meta.data$celltype)
names(pheno) <- colnames(seu2)
plotCytoTRACE(results, phenotype = pheno)
# Get umap embeddings from Seurat...
# plotCytoTRACE(results, phenotype = pheno, gene = "TIMP1", emb = umapEmb)
