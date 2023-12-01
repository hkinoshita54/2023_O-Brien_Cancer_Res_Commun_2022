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
options(Seurat.object.assay.version = "v5")
library(monocle3)

seu2 <- readRDS(file = "RDSfiles/seu2.RDS")

# re-cluster cells from 12weeks and freeze-thaw within seu2 ----
seu12wft <- subset(seu2, subset = orig.ident %in% c("1", "2", "3", "4"))
seu12wft <- NormalizeData(seu12wft, verbose = FALSE)
seu12wft <- FindVariableFeatures(seu12wft, verbose = FALSE)
seu12wft <- ScaleData(seu12wft, verbose = FALSE)
seu12wft <- RunPCA(seu12wft, npcs = 10, verbose = FALSE)
seu12wft <- FindNeighbors(seu12wft, dims = 1:10, verbose = FALSE)
seu12wft <- FindClusters(seu12wft, resolution = 0.5, verbose = FALSE)
seu12wft <- RunUMAP(seu12wft, dims = 1:10, verbose = FALSE)
DimPlot(seu12wft, label = TRUE, repel = TRUE) + NoAxes()
DimPlot(seu12wft, group.by = "orig.ident") + NoAxes()
DimPlot(seu12wft, group.by = "Kras_status") + NoAxes()
DimPlot(seu12wft, group.by = "infection") + NoAxes()
DimPlot(seu12wft, split.by = "orig.ident", ncol = 2) + NoAxes()

FeaturePlot(seu12wft,features = "Epcam", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu12wft,features = "Ptprc", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu12wft,features = "Dcn", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu12wft,features = "Muc5ac", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu12wft,features = "Muc6", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu12wft,features = "Muc4", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu12wft,features = "Gif", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu12wft,features = "Atp4b", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu12wft,features = "Mki67", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu12wft,features = "Ptgs2", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()

seu12wft <- RenameIdents(seu12wft, 
                         `0` = "Pit3", `1` = "Pit1", `2` = "Parietal1", 
                         `3` = "Neck2", `4` = "Pit2", `5` = "Ki67", 
                         `6` = "Chief", `7` = "Neck1", `8` = "Parietal2", `9` = "Immune")
seu12wft$celltype <- Idents(seu12wft)
DimPlot(seu12wft, label = TRUE, repel = TRUE) + NoAxes()

# selecting cells by using monocle3 function choose_cells()
cds <- as.cell_data_set(seu12wft)
cds <- choose_cells(cds)
seu12wft <- seu12wft[,cds@colData@rownames]
seu12wft <- subset(seu12wft, subset = celltype %in% c("Parietal1", "Parietal2", "Immune"), invert = TRUE)
DimPlot(seu12wft, label = TRUE, repel = TRUE) + NoAxes()
DimPlot(seu12wft, split.by = "orig.ident", ncol = 2) + NoAxes()

FeaturePlot(seu12wft,features = "Muc5ac", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu12wft,features = "Muc6", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu12wft,features = "Muc4", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu12wft,features = "Gif", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu12wft,features = "Atp4b", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu12wft,features = "Mki67", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu12wft,features = "Stmn1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu12wft,features = "Cd44", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()

# save data for RNA velocity ----
# save metadata table
seu12wft$barcode <- colnames(seu12wft)
seu12wft$UMAP_1 <- seu12wft@reductions$umap@cell.embeddings[,1]
seu12wft$UMAP_2 <- seu12wft@reductions$umap@cell.embeddings[,2]
write.csv(seu12wft@meta.data, file='out/seu12wft_metadata.csv', quote=F, row.names=F)

# write expression counts matrix
library(Matrix)
counts_matrix <- LayerData(seu12wft, assay='RNA', layer='counts')
writeMM(counts_matrix, file='out/seu12wft_counts.mtx')

# write dimesnionality reduction matrix, in this example case pca matrix
write.csv(seu12wft@reductions$pca@cell.embeddings, file='out/seu12wft_pca.csv', quote=F, row.names=F)

# write gene names
write.table(
  data.frame('gene'=rownames(counts_matrix)),file='out/seu12wft_gene_names.csv',
  quote=F,row.names=F,col.names=F
)

saveRDS(seu12wft, file = "RDSfiles/seu12wft.RDS")

# CytoTRACE ----
library(CytoTRACE)
counts_matrix <- LayerData(seu12wft, assay='RNA', layer='counts') %>% as.data.frame()
obj_cell_type_anno <- as.data.frame(seu12wft@meta.data$celltype)
results <- CytoTRACE(counts_matrix, ncores = 4)
pheno <- as.character(seu12wft@meta.data$celltype)
names(pheno) <- colnames(seu12wft)
plotCytoTRACE(results, phenotype = pheno)

# monocle3 ----
# https://rpubs.com/mahima_bose/Seurat_and_Monocle3_p
cds <- as.cell_data_set(seu12wft)
fData(cds)$gene_short_name <- rownames(fData(cds))

# assign partitions
recreate.partitions <- c(rep(1, length(cds@colData@rownames)))
names(recreate.partitions) <- cds@colData@rownames
recreate.partitions <- as.factor(recreate.partitions)
cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partitions

# assign clusters
list.cluster <- seu12wft@active.ident
cds@clusters@listData[["UMAP"]][["clusters"]] <- list.cluster

# assign UMAP coordinates
cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- seu12wft@reductions$umap@cell.embeddings

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

# find genes that change as a function of pseudotime
# deg <- graph_test(cds, neighbor_graph = "principal_graph")
# deg %>% arrange(q_value) %>% filter(status == "OK") %>% head()
# FeaturePlot(seu12wft, features = c("Cox5b", "Zap70", "Ptma", "Rgs1"))

# add pseudotime values to the seurat object
# seu12wft$pseudotime <- pseudotime(cds)
# FeaturePlot(seu12wft, features = "pseudotime") + NoAxes()
# RidgePlot(seu12wft, features = c("pseudotime"), sort = TRUE, idents = c("Ki67", "Neck1", "Neck2", "Chief"))
# RidgePlot(seu12wft, features = c("pseudotime"), sort = TRUE, idents = c("Ki67", "Pit1", "Pit2", "Pit3"))
# my_genes <- row.names(subset(fData(cds), gene_short_name %in% c("Mki67", "Muc6", "Pgc", "Muc4", "Muc5ac"))) 
# cds_subset <- cds[my_genes,]
# plot_genes_in_pseudotime(cds_subset, color_cells_by = "monocle3_pseudotime" )



