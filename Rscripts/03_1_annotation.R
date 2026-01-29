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
options(Seurat.object.assay.version = "v5")


####################
# load data ----
seu <- readRDS(file = "RDSfiles/seu_05_lognorm_harmony.RDS")

####################
# annotation ----
seu <- RenameIdents(
  seu, 
  `0` = "Pit1", `1` = "Pit2", `2` = "Tcell_Nkg7", `3` = "BEC", `4` = "Pit3", 
  `5` = "Macrophage_Ccr2", `6` = "Macrophage_C1qb", `7` = "LEC", `8` = "Tuft", `9` = "LEC", 
  `10` = "Squamous", `11` = "Tcell_Rag1", `12` ="cDC", `13` = "Parietal", `14` = "Neck1", 
  `15` = "Neck2", `16` = "Granulocyte", `17` = "Squamous", `18` = "Stem", `19` = "Stem", `20` = "Endocrine1",
  `21` = "Endocrine2", `22` = "Bcell", `23` = "Erythrocyte", `24` = "Myocyte", `25` = "Serosa",
  `26` = "Tuft", `27` = "Fibroblast", `28` = "BEC", `29` = "Chief", `30` = "pDC", `31` = "Tcell_Rag1"
)
seu$celltype <- Idents(seu)
seu$celltype <- factor(seu$celltype, levels = c("Stem", "Pit1", "Pit2", "Pit3", "Neck1", "Neck2", "Chief", "Parietal", "Endocrine1", "Endocrine2","Tuft", "Squamous",
                                                "Tcell_Rag1", "Tcell_Nkg7", "Bcell", "Macrophage_C1qb", "Macrophage_Ccr2", "Granulocyte", "cDC", "pDC", "Erythrocyte",
                                                "BEC", "LEC", "Fibroblast", "Myocyte", "Serosa"))
Idents(seu) <- seu$celltype
DimPlot(seu, label = TRUE, repel = TRUE) + NoAxes()
ggsave("annotation_manual.png", path = "plots/04_lognorm_harmony/annotation", width = 8, height = 6, units = "in", dpi = 150)


####################
# recheck feature plot ----
files <- list.files(path = "gene_set/annotation/", full.names = TRUE)
features <- lapply(files, FUN = readLines) %>% unlist()
for(i in 1:length(features)){
  p <- FeaturePlot(seu, features = features[i], cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
  ggsave(paste0(as.character(features[i]), ".png"), plot = p, path = "plots/04_lognorm_harmony/annotation", width = 5, height = 5, units = "in", dpi = 150)
}


saveRDS(seu, file = "RDSfiles/seu_05_lognorm_harmony.RDS")


####################
# check dot plot ----
epi_features = readLines("gene_set/annotation/01_epi_markers.txt")
epi <- subset(seu, idents = c("Stem", "Pit1", "Pit2", "Pit3", "Neck1", "Neck2", "Chief", "Parietal", "Endocrine1", "Endocrine2","Tuft", "Squamous"))
DotPlot(epi, features = epi_features) + RotatedAxis()
# ggsave("dotplot_epi.png", path = "plots/04_lognorm_harmony/annotation", width = 8, height = 6, units = "in", dpi = 150)

imm_features = readLines("gene_set/annotation/02_imm_markers.txt")
imm <- subset(seu, idents = c("Tcell_Rag1", "Tcell_Nkg7", "Bcell", "Macrophage_C1qb", "Macrophage_Ccr2", "Granulocyte", "cDC", "pDC"))
DotPlot(imm, features = imm_features) + RotatedAxis()
# ggsave("dotplot_imm.png", path = "plots/04_lognorm_harmony/annotation", width = 8, height = 6, units = "in", dpi = 150)

str_features = readLines("gene_set/annotation/03_str_markers.txt")
str <- subset(seu, idents = c("BEC", "LEC", "Fibroblast", "Myocyte", "Serosa"))
DotPlot(str, features = str_features) + RotatedAxis()
# ggsave("dotplot_str.png", path = "plots/04_lognorm_harmony/annotation", width = 8, height = 6, units = "in", dpi = 150)