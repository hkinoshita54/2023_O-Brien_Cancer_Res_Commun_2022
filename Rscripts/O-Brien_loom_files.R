# load packages ----
library(tidyverse)
library(cowplot)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(Seurat)
library(SeuratWrappers)
library(SeuratDisk)
options(Seurat.object.assay.version = "v3") # needed to run as.Seurat()

# load loom files as seurat object, merge all the samples, then convert it to anndata ----
# https://htmlpreview.github.io/?https://github.com/satijalab/seurat-wrappers/blob/master/docs/scvelo.html
sra_table <- read.delim(file = "data/SraRunTable.txt", sep = ",") # read SraRunTable from GSE224840
ids <- as.character(sra_table$Sample.Name) # get character vector of sample ids

ldat_list <- list()
for(id in ids){
  # load loom file, add sample id to the seurat object
  file <- dir(path = "./data/loom/", pattern = paste0("^", id, "_"), full.names = TRUE)
  ldat <- ReadVelocity(file = file)
  ldat <- as.Seurat(x = ldat)
  ldat$orig.ident <- id
  
  # change cell names according to cellranger output
  colnames(ldat) <- sapply(strsplit(colnames(ldat), ":"), "[", 2)
  colnames(ldat) <- gsub("x", "-1", colnames(ldat))
  # colnames(ldat) <- paste(id, colnames(ldat), sep = "_")
  
  # append ldat to ldat_list
  ldat_list <- append(ldat_list, ldat)
}
ldat <- merge(x = ldat_list[[1]], y = ldat_list[2:length(ids)], add.cell.ids = ids)

saveRDS(ldat, file = "RDSfiles/ldat.RDS")

# generate anndata from loom file, only selecting cells in the processed seurat object ----
seu2 <- readRDS(file = "RDSfiles/seu2.RDS")
ldat2 <- ldat[, colnames(seu2)]
ldat2[["RNA"]] <- ldat2[["spliced"]]
DefaultAssay(ldat2) <- "RNA"
SaveH5Seurat(ldat2, filename = "out/ldat2.h5Seurat", overwrite = TRUE)
Convert("out/ldat2.h5Seurat", dest = "h5ad", overwrite = TRUE)

# generate anndata from loom file, selecting cells and features in the processed seurat object ----
seu2 <- readRDS(file = "RDSfiles/seu2.RDS")
ldat3 <- ldat[rownames(seu2), colnames(seu2)]
ldat3[["RNA"]] <- ldat3[["spliced"]]
DefaultAssay(ldat3) <- "RNA"
SaveH5Seurat(ldat3, filename = "out/ldat3.h5Seurat", overwrite = TRUE)
Convert("out/ldat3.h5Seurat", dest = "h5ad", overwrite = TRUE)

# generate anndata from loom file, selecting cells and features in the processed seurat object ----
ldat3 <- ldat[rownames(seu3_2), colnames(seu3_2)]
ldat3[["RNA"]] <- ldat3[["spliced"]]
DefaultAssay(ldat3) <- "RNA"
SaveH5Seurat(ldat3, filename = "out/ldat3_2.h5Seurat", overwrite = TRUE)
Convert("out/ldat3_2.h5Seurat", dest = "h5ad", overwrite = TRUE)

# generate anndata from loom file, selecting cells in the processed seurat object of 12 weeks and freeze/thaw ----
seu12wft <- readRDS(file = "RDSfiles/seu12wft.RDS")
ldat12wft <- ldat[rownames(seu2), colnames(seu2)]
ldat12wft[["RNA"]] <- ldat12wft[["spliced"]]
DefaultAssay(ldat12wft) <- "RNA"
SaveH5Seurat(ldat12wft, filename = "out/ldat12wft.h5Seurat", overwrite = TRUE)
Convert("out/ldat12wft.h5Seurat", dest = "h5ad", overwrite = TRUE)