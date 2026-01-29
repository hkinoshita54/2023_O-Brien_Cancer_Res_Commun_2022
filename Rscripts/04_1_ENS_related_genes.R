####################
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


####################
# load data and check features by Vln ----
seu <- readRDS(file = "RDSfiles/seu_05_lognorm_harmony.RDS")


####################
# check genes of interest from Dr. Hayakawa ----
features = readLines("gene_set/ENS_related.txt")
for (i in 1:length(features)){
  tryCatch({
    p <- FeaturePlot(seu, features = features[i], cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE)
    + NoAxes() + NoLegend()
    ggsave(paste0(as.character(features[i]), ".png"), plot = p, 
           path = "plots/04_lognorm_harmony/annotation", 
           width = 5, height = 5, units = "in", dpi = 150)
  }, error = function(e){cat("ERROR :", conditionMessage(e), "\n")})
}

epi <- subset(seu, idents = c("Stem", "Pit1", "Pit2", "Pit3", "Neck1", "Neck2", "Chief", "Parietal", "Endocrine1", "Endocrine2","Tuft", "Squamous"))
DotPlot(epi, features = features) + RotatedAxis()

imm <- subset(seu, idents = c("Tcell_Rag1", "Tcell_Nkg7", "Bcell", "Macrophage_C1qb", "Macrophage_Ccr2", "Granulocyte", "cDC", "pDC"))
DotPlot(imm, features = features) + RotatedAxis()

str <- subset(seu, idents = c("BEC", "LEC", "Fibroblast", "Myocyte", "Serosa"))
DotPlot(str, features = features) + RotatedAxis()
