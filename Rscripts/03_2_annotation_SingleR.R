####################
# load packages ----
library(tidyverse)
library(cowplot)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(Seurat)
library(SeuratWrappers)
library(SeuratDisk)
# options(Seurat.object.assay.version = "v5")



####################
# load data ----
seu <- readRDS(file = "RDSfiles/seu_05_lognorm_harmony.RDS")


####################
# automated annotation by SingleR and celldex ----
library(SingleR)
library(celldex)
options(Seurat.object.assay.version = "v3")
seu[["RNA"]] <- as(object = seu[["RNA"]], Class = "Assay")
sce <- as.SingleCellExperiment(seu)
# ref <- MouseRNAseqData()
ref <- ImmGenData()
res <- SingleR(test = sce, ref = ref, assay.type.test=1, labels = ref$label.main)
seu$SingleR.labels <- res$labels
DimPlot(seu, group.by = "SingleR.labels") + NoAxes()


