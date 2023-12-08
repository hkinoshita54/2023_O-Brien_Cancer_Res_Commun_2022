library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
options(Seurat.object.assay.version = "v3") # needed to run as.Seurat()

curl::curl_download(url = 'http://pklab.med.harvard.edu/velocyto/mouseBM/SCG71.loom', destfile　= '~/Downloads/SCG71.loom')
ldat <- ReadVelocity(file = "~/Downloads/SCG71.loom")
bm <- as.Seurat(x = ldat)
bm[["RNA"]] <- bm[["spliced"]]
bm <- SCTransform(bm)
bm <- RunPCA(bm)
bm <- RunUMAP(bm, dims = 1:20)
bm <- FindNeighbors(bm, dims = 1:20)
bm <- FindClusters(bm)
DefaultAssay(bm) <- "RNA"
SaveH5Seurat(bm, filename = "mouseBM.h5Seurat")
Convert("mouseBM.h5Seurat", dest = "h5ad")