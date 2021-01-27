library(Seurat)
library(tidyverse)
set.seed(2021)

sc <- readRDS("../data/seurat.normalized.rds")

sc <- RunPCA(sc, npcs=100, features = VariableFeatures(object = sc)[1:1000])
sc <- FindNeighbors(sc, dims=1:50, nn.eps=0.1)
sc <- FindClusters(sc, resolution=0.15, n.start=10)
sc <- RunUMAP(sc, dims=1:50, min.dist=0.15, n.neighbors=30)

saveRDS(sc, file="../data/seurat.processed.rds")