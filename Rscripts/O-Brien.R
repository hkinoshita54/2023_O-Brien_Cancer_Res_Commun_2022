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


# load data ----
sra_table <- read.delim(file = "data/SraRunTable.txt", sep = ",") # read SraRunTable from GSE224840
ids <- as.character(sra_table$Sample.Name) # get character vector of sample ids

seu_list <- list()
for (id in ids){
  mtx <- dir(path = "./data/GSE224840_RAW/", pattern = paste0("matrix", id, "\\.mtx"), full.names = TRUE)
  features <- dir(path = "./data/GSE224840_RAW/", pattern = paste0("features", id, "\\.tsv"), full.names = TRUE)
  barcodes <- dir(path = "./data/GSE224840_RAW/", pattern = paste0("barcodes", id, "\\.tsv"), full.names = TRUE)
  cts <- ReadMtx(mtx = mtx, features = features, cells = barcodes)
  seu <- CreateSeuratObject(counts = cts, project = id, min.cells = 3, min.features = 200)
  seu_list <- append(seu_list, seu)
}
seu <- merge(x = seu_list[[1]], y = seu_list[2:length(ids)], add.cell.ids = ids)

# meta_data <- read.delim("data/GSE224840_scRNA_fullset_metadata.csv.gz", sep = ",", stringsAsFactors = F, header = T, row.names = 1)
# meta_data <- select(meta_data, Infection, Kras.Status., Time..weeks., Preparation, label, SPEM)
# colnames(meta_data) <- c("Infection", "Kras_tatus", "time_weeks", "Preparation", "label", "SPEM")
####
#### cells names in meta_data not match with cell names in seurat object
#### couldn't solve the problem, proceed without meta_data

seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^mt-")
# VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
seu <- subset(seu, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 25)

# cluster without integration ----
seu <- JoinLayers(seu)
seu <- NormalizeData(seu, verbose = FALSE)
seu <- FindVariableFeatures(seu, verbose = FALSE)
seu <- ScaleData(seu, verbose = FALSE)

# seu$orig.ident <- factor(seu$orig.ident)
# seu <- SCTransform(seu, vars.to.regress = c("percent.mt", "orig.ident"))

seu <- RunPCA(seu, npcs = 20, verbose = FALSE)
seu <- FindNeighbors(seu, dims = 1:20, verbose = FALSE)
seu <- FindClusters(seu, resolution = 1.5, verbose = FALSE)
seu <- RunUMAP(seu, dims = 1:20, verbose = FALSE)
DimPlot(seu, label = TRUE, repel = TRUE) + NoAxes()
DimPlot(seu, group.by = "orig.ident") + NoAxes()

FeaturePlot(seu,features = "Epcam", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "Ptprc", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "Dcn", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "Muc5ac", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "Muc6", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "Muc4", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "Gif", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "Atp4b", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "Tff3", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "Chga", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "Dclk1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "Krt5", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "Hba-a1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()

saveRDS(seu, file = "RDSfiles/seu.RDS")

# cluster with integration ----
# seu <- NormalizeData(seu, verbose = FALSE)
# seu <- FindVariableFeatures(seu, verbose = FALSE)
# seu <- ScaleData(seu, verbose = FALSE)
# seu <- RunPCA(seu, npcs = 30, verbose = FALSE)
# 
# seu <- IntegrateLayers(
#   object = seu, method = scVIIntegration,
#   new.reduction = "integrated.scvi",
#   conda_env = "/Users/kinoshitahiroto/miniconda3/envs/scvi-env",
#   verbose = TRUE
# )
# 
# seu <- FindNeighbors(seu, reduction = "integrated.scvi", dims = 1:30, verbose = FALSE)
# seu <- FindClusters(seu, resolution = 2, cluster.name = "scvi_clusters", verbose = FALSE)
# seu <- RunUMAP(seu, reduction = "integrated.scvi", dims = 1:30, reduction.name = "umap.scvi", verbose = FALSE)
# 
# DimPlot(seu, label = TRUE, repel = TRUE) + NoAxes()
# DimPlot(seu, label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
# DimPlot(seu, group.by = "orig.ident") + NoAxes()


# subset pit, neck, chief, parietal ----
## lognormalize without integration ----
seu2 <- subset(seu, idents = c(0,3,6,7,14,15,19,20,28,29))
seu2 <- NormalizeData(seu2, verbose = FALSE)
seu2 <- FindVariableFeatures(seu2, verbose = FALSE)
seu2 <- ScaleData(seu2, verbose = FALSE)
seu2 <- RunPCA(seu2, npcs = 10, verbose = FALSE)
seu2 <- FindNeighbors(seu2, dims = 1:10, verbose = FALSE)
seu2 <- FindClusters(seu2, resolution = 1, verbose = FALSE)
seu2 <- RunUMAP(seu2, dims = 1:10, verbose = FALSE)
DimPlot(seu2, label = TRUE, repel = TRUE) + NoAxes()
DimPlot(seu2, group.by = "orig.ident") + NoAxes()

FeaturePlot(seu2,features = "Epcam", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu2,features = "Ptprc", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu2,features = "Dcn", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu2,features = "Muc5ac", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu2,features = "Muc6", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu2,features = "Muc4", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu2,features = "Gif", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu2,features = "Atp4b", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu2,features = "Mki67", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu2,features = "Clu", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu2,features = "Ly6a", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu2,features = "Hba-a1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()

Idents(seu2) <- "seurat_clusters"
seu2 <- RenameIdents(seu2, 
                           `0` = "Pit0", `1` = "Pit1", `2` = "Pit2",`3` = "Pit3", `4` = "Pit4", `5` = "Neck1", 
                           `6` = "Neck2", `7` = "Pit5", `8` = "Ki67", `9` = "Uncomitted", `10` = "Parietal1", `11` = "Parietal2", `12` ="Pit6", `13` = "Chief", `14` = "Parietal3", `15` = "Pit7", `16` = "??")
seu2$celltype <- Idents(seu2)
seu2$celltype <- factor(seu2$celltype, levels = c("Ki67","Pit0","Pit1","Pit2","Pit3","Pit4","Pit5","Pit6","Pit7","Neck1","Neck2","Chief","Parietal1","Parietal2","Parietal3","Uncomitted","??"))
Idents(seu2) <- seu2$celltype
DimPlot(seu2, label = TRUE, repel = TRUE) + NoAxes()

# add some more meta data to seu2 ----

# Kras status
Idents(seu2) <- "orig.ident"
Kras <- WhichCells(seu2, ident = c("3","4","2204","2266","2547","2580"))
WT <- WhichCells(seu2, ident = c("1","2"))
seu2 <- SetIdent(seu2, cells = Kras, value = "mut")
seu2$Kras_status <- Idents(seu2)
seu2 <- SetIdent(seu2, cells = WT, value = "wt")
seu2$Kras_status <- Idents(seu2)

Idents(seu2) <- "orig.ident"
Hp <- WhichCells(seu2, ident = c("2","4","2204","2547"))
mock <- WhichCells(seu2, ident = c("1","3","2266","2580"))
seu2 <- SetIdent(seu2, cells = Hp, value = "Hp")
seu2$infection <- Idents(seu2)
seu2 <- SetIdent(seu2, cells = mock, value = "mock")
seu2$infection <- Idents(seu2)

Idents(seu2) <- "orig.ident"
fresh <- WhichCells(seu2, ident = c("2547","2580"))
fr_th <- WhichCells(seu2, ident = c("1","2","3","4","2204","2266"))
seu2 <- SetIdent(seu2, cells = fresh, value = "fresh")
seu2$preparation <- Idents(seu2)
seu2 <- SetIdent(seu2, cells = fr_th, value = "freeze_thaw")
seu2$preparation <- Idents(seu2)

Idents(seu2) <- "orig.ident"
b1 <- WhichCells(seu2, ident = c("2266","2204"))
b2 <- WhichCells(seu2, ident = c("2580","2547"))
b3 <- WhichCells(seu2, ident = c("1","2","3","4"))
seu2 <- SetIdent(seu2, cells = b1, value = "batch1")
seu2$gem_prep <- Idents(seu2)
seu2 <- SetIdent(seu2, cells = b2, value = "batch2")
seu2$gem_prep <- Idents(seu2)
seu2 <- SetIdent(seu2, cells = b3, value = "batch3")
seu2$gem_prep <- Idents(seu2)

Idents(seu2) <- "celltype"
saveRDS(seu2, file = "RDSfiles/seu2.RDS")

# SaveH5Seurat(seu2, filename = "out/seu2.h5Seurat", overwrite = TRUE)
# Convert("out/seu2.h5Seurat", dest = "h5ad", overwrite = TRUE)
### Convert() did not work, maybe d/t version of Seurat or v5 objects, although it worked in the other Rscript "O-Brien_velocity_test3.R"

## convert Seurat object to anndata manually following the tutorial below ----
# https://smorabit.github.io/tutorials/8_velocyto/

# save metadata table
seu2$barcode <- colnames(seu2)
seu2$UMAP_1 <- seu2@reductions$umap@cell.embeddings[,1]
seu2$UMAP_2 <- seu2@reductions$umap@cell.embeddings[,2]
write.csv(seu2@meta.data, file='out/seu2_metadata.csv', quote=F, row.names=F)

# write expression counts matrix
library(Matrix)
counts_matrix <- LayerData(seu2, assay='RNA', layer='counts')
writeMM(counts_matrix, file='out/seu2_counts.mtx')

# write dimesnionality reduction matrix, in this example case pca matrix
write.csv(seu2@reductions$pca@cell.embeddings, file='out/seu2_pca.csv', quote=F, row.names=F)

# write gene names
write.table(
  data.frame('gene'=rownames(counts_matrix)),file='out/seu2_gene_names.csv',
  quote=F,row.names=F,col.names=F
)