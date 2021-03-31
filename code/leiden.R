library(Seurat)

# options(future.globals.maxSize=4e5*1024^2)
# sc <- readRDS("../data/seurat.normalized.rds")
# sc <- RunPCA(sc, npcs=100, features = VariableFeatures(object = sc)[1:1000])
# sc <- FindNeighbors(sc, dims=1:50, nn.eps=0.1)
sc <- readRDS("../data/seurat.neighbors.rds")
sc <- FindClusters(sc, resolution=0.35, n.start=10)
saveRDS(sc, "../data/seurat.louvain055.rds")

sc.sub <- sc[,!sc$seurat_clusters %in% c(11, 12)]
sc.sub <- sc.sub %>%
  ScaleData() %>%
  RunPCA(npcs=100, features=VariableFeatures(object=sc.sub)[1:1000])
saveRDS(sc.sub, "../data/seurat.sub.leiden.rds")
  