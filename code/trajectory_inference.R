library(tidyverse)
library(Seurat)
library(slingshot)
library(brewer.pal)

sc_cm <- readRDS("../data/seurat.cm.rds")

ti_reduced <- sc_cm@reductions$pca@cell.embeddings[,1:3]

slung <- slingshot(ti_reduced, clusterLabels=as.character(sc_cm$type), 
                   start.clus="IPSC", end.clus="CM",
                   approx_points=100)
sc_cm$pseudotime <- slingPseudotime(slung)

saveRDS(sc_cm, "../data/seurat.cm.rds")
saveRDS(slung, "../data/slingshot.cm.rds")


sc_cf <- readRDS("../data/seurat.cf.rds")

cf_reduced <- sc_cf@reductions$pca@cell.embeddings[,1:3]

slung <- slingshot(cf_reduced, clusterLabels=as.character(sc_cf$type), 
                   start.clus="IPSC", end.clus="CF",
                   approx_points=100)
sc_cf$pseudotime <- slingPseudotime(slung)

saveRDS(sc_cm, "../data/seurat.cm.rds")
saveRDS(slung, "../data/slingshot.cm.rds")