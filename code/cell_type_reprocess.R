library(Seurat)
library(tidyverse)
set.seed(2021)

sc <- readRDS("../data/seurat.normalized.rds")

conservative_omit <- read_tsv("../data/outliers_conservative.tsv")$cell
sc_cons <- sc[,setdiff(colnames(sc), conservative_omit)]
sc_cons <- ScaleData(sc_cons)
sc_cons <- RunPCA(sc_cons, npcs=100, features = VariableFeatures(object = sc_cons)[1:1000])
sc_cons <- FindNeighbors(sc_cons, dims=1:50)
sc_cons <- FindClusters(sc_cons, resolution=0.15, n.start=10)
# sc_cons$lo_res_clusters <- sc_cons$seurat_clusters
# sc_cons <- FindClusters(sc_cons, resolution=0.65, n.start=10)
# sc_cons$hi_res_clusters <- sc_cons$seurat_clusters
sc_cons <- RunUMAP(sc_cons, dims=1:50, min.dist=0.15, n.neighbors=100)

liberal_omit <- read_tsv("../data/outliers_liberal.tsv")$cell
sc_lib <- sc[,setdiff(colnames(sc), liberal_omit)]
sc_lib <- ScaleData(sc_lib)
sc_lib <- RunPCA(sc_lib, npcs=100, features = VariableFeatures(object = sc_lib)[1:1000])
sc_lib <- FindNeighbors(sc_lib, dims=1:50, nn.eps=0.1)
sc_lib <- FindClusters(sc_lib, resolution=0.15, n.start=10)
# sc_lib$lo_res_clusters <- sc_lib$seurat_clusters
# sc_lib <- FindClusters(sc_lib, resolution=0.65, n.start=10)
# sc_lib$hi_res_clusters <- sc_lib$seurat_clusters
sc_lib <- RunUMAP(sc_lib, dims=1:50, min.dist=0.15, n.neighbors=50)

saveRDS(sc_cons, file="../data/seurat.conservative.rds")
saveRDS(sc_lib, file="../data/seurat.liberal.rds")