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
sra_table <- read.delim(file = "data/SraRunTable.txt", sep = ",") # read SraRunTable from GSE224840
ids <- sra_table$Sample.Name %>% sort() %>% as.character() # get character vector of sample ids

seu_list <- list()
for (i in 1:length(ids)){
  mtx <- dir(path = "./data/GSE224840_RAW/", pattern = paste0("matrix", ids[i], "\\.mtx"), full.names = TRUE)
  features <- dir(path = "./data/GSE224840_RAW/", pattern = paste0("features", ids[i], "\\.tsv"), full.names = TRUE)
  barcodes <- dir(path = "./data/GSE224840_RAW/", pattern = paste0("barcodes", ids[i], "\\.tsv"), full.names = TRUE)
  cts <- ReadMtx(mtx = mtx, features = features, cells = barcodes)
  seu <- CreateSeuratObject(counts = cts, project = ids[i], min.cells = 3, min.features = 200)
  seu_list[i] <- seu
}
seu <- merge(x = seu_list[[1]], y = seu_list[2:length(ids)], add.cell.ids = ids)


####################
#### QC: filter by percent.mt and nFeature_RNA ####
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^mt-")
VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
seu <- subset(seu, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 25)

####################
#### add meta data ####

# meta data from supplementary file of GEO
meta_data <- read.delim("data/GSE224840_scRNA_fullset_metadata.csv.gz", sep = ",", stringsAsFactors = F, header = T, row.names = 1)
all(rownames(meta_data) %in% colnames(seu))
#### cell names in meta_data do not match cell names in seurat object, couldn't solve this problem, proceed without meta_data

# change levels of orig.ident
seu$orig.ident <- factor(seu$orig.ident, levels = c("1", "2", "3", "4", "2204", "2266", "2547", "2580"))

# Kras status
Idents(seu) <- "orig.ident"
Kras <- WhichCells(seu, ident = c("3","4","2204","2266","2547","2580"))
WT <- WhichCells(seu, ident = c("1","2"))
seu <- SetIdent(seu, cells = Kras, value = "mut")
seu$Kras_status <- Idents(seu)
seu <- SetIdent(seu, cells = WT, value = "wt")
seu$Kras_status <- Idents(seu)

# Hp infection
Idents(seu) <- "orig.ident"
Hp <- WhichCells(seu, ident = c("2","4","2204","2547"))
mock <- WhichCells(seu, ident = c("1","3","2266","2580"))
seu <- SetIdent(seu, cells = Hp, value = "Hp")
seu$infection <- Idents(seu)
seu <- SetIdent(seu, cells = mock, value = "mock")
seu$infection <- Idents(seu)

# preparation (fresh or freeze/thaw)
Idents(seu) <- "orig.ident"
fresh <- WhichCells(seu, ident = c("2547","2580"))
fr_th <- WhichCells(seu, ident = c("1","2","3","4","2204","2266"))
seu <- SetIdent(seu, cells = fresh, value = "fresh")
seu$preparation <- Idents(seu)
seu <- SetIdent(seu, cells = fr_th, value = "freeze_thaw")
seu$preparation <- Idents(seu)

# Batch (GEM preparation)
Idents(seu) <- "orig.ident"
b1 <- WhichCells(seu, ident = c("2266","2204"))
b2 <- WhichCells(seu, ident = c("2580","2547"))
b3 <- WhichCells(seu, ident = c("1","2","3","4"))
seu <- SetIdent(seu, cells = b1, value = "batch1")
seu$gem_prep <- Idents(seu)
seu <- SetIdent(seu, cells = b2, value = "batch2")
seu$gem_prep <- Idents(seu)
seu <- SetIdent(seu, cells = b3, value = "batch3")
seu$gem_prep <- Idents(seu)


saveRDS(seu, file = "RDSfiles/seu_01_mt25.RDS")