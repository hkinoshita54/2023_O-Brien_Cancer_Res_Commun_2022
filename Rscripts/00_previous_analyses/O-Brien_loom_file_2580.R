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
id = "2580"
file <- dir(path = "./data/loom/", pattern = paste0("^", id, "_"), full.names = TRUE)
ldat <- ReadVelocity(file = file)
ldat <- as.Seurat(x = ldat)
ldat$orig.ident <- id

# change cell names according to cellranger output
colnames(ldat) <- sapply(strsplit(colnames(ldat), ":"), "[", 2)
colnames(ldat) <- gsub("x", "-1", colnames(ldat))
colnames(ldat) <- paste(id, colnames(ldat), sep = "_")

# generate anndata from loom file, only selecting cells in the processed seurat object ----
seu2580 <- readRDS(file = "RDSfiles/seu2580.RDS")
ldat2580 <- ldat[, colnames(seu2580)]
ldat2580[["RNA"]] <- ldat2580[["spliced"]]
DefaultAssay(ldat2580) <- "RNA"
SaveH5Seurat(ldat2580, filename = "out/ldat2580.h5Seurat", overwrite = TRUE)
Convert("out/ldat2580.h5Seurat", dest = "h5ad", overwrite = TRUE)