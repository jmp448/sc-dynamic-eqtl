library(Seurat)
library(tidyverse)
library(Matrix)
library(reticulate)
use_condaenv("singlecell")

sc <- readRDS("../data/seurat.pca.rds")

# cluster genes
feats <- sc@reductions$pca@feature.loadings
feats_snn <- FindNeighbors(feats)
feats_clusts <- as_tibble(FindClusters(feats_snn$snn, resolution=0.5, algorithm=4), rownames="gene") %>%
  `colnames<-`(c("gene", "cluster"))
feats_umap <- RunUMAP(feats_snn$snn)@cell.embeddings

# manually separate the cell cycle genes from some mesoderm stage genes 
feats_clusts$cluster <- as.integer(feats_clusts$cluster)
feats_clusts[(feats_umap[,1] < -20) & (feats_clusts$cluster==5), "cluster"] = 9
feats_clusts$cluster <- factor(feats_clusts$cluster)

feats_embedding <- as_tibble(feats_umap, rownames="gene") %>%
  inner_join(feats_clusts, by="gene")

ggplot(feats_embedding, aes(x=UMAP_1, y=UMAP_2, color=cluster)) +
  geom_point() +
  theme_classic()

# recompute pcs using just the interesting genes
differentiation.set <- filter(feats_clusts, cluster %in% c(1,2,5,7))$gene

sc.diff <- RunPCA(sc, features=differentiation.set)
PCAPlot(sc.diff, group.by="diffday")

# get module-specific scores 
c1 <- filter(feats_clusts, cluster==1)$gene
c2 <- filter(feats_clusts, cluster==2)$gene
c5 <- filter(feats_clusts, cluster==5)$gene
c7 <- filter(feats_clusts, cluster==7)$gene
c9 <- filter(feats_clusts, cluster==9)$gene

c1score <- RunPCA(sc, features=c1)@reductions$pca@`cell.embeddings`[,1] 
sc$c1_pcscore <- c1score

c2score <- RunPCA(sc, features=c2)@reductions$pca@`cell.embeddings`[,1] 
sc$c2_pcscore <- c2score

c5score <- RunPCA(sc, features=c5)@reductions$pca@`cell.embeddings`[,1] 
sc$c5_pcscore <- c5score

c7score <- RunPCA(sc, features=c7) 
sc$c7_pcscore <- c7score@reductions$pca@`cell.embeddings`[,1]

c9score <- RunPCA(sc, features=c9) 
sc$c9_pcscore <- c9score@reductions$pca@`cell.embeddings`[,1]

# look at each module score vs real differentiation time
sc$day <- as.numeric(sc$diffday)
FeatureScatter(sc, "day", "c1_pcscore", group.by="diffday") #cm
FeatureScatter(sc, "day", "c2_pcscore", group.by="diffday") #cf
FeatureScatter(sc, "day", "c5_pcscore", group.by="diffday") #mesoderm
FeatureScatter(sc, "day", "c7_pcscore", group.by="diffday") #mesoderm
FeatureScatter(sc, "day", "c9_pcscore") #iPSC/ cell cycle
saveRDS(sc, "../data/seurat.genemodules.rds")

# load the full embedded data
sc.embedded <- readRDS("../data/seurat.annotated.rds")
cell.scores <- as_tibble(sc@meta.data, rownames="cell") %>% select(cell, ends_with("pcscore"))
cell.metadata <- tibble("cell"=colnames(sc.embedded)) %>% left_join(cell.scores, by="cell") %>% column_to_rownames("cell")
sc.embedded <- AddMetaData(sc.embedded, cell.metadata)
sc.embedded$c2_pcscore = -1 * sc.embedded$c2_pcscore
sc.embedded$c5_pcscore = -1 * sc.embedded$c5_pcscore
DimPlot(sc.embedded, reduction="fa", group.by="type")
FeaturePlot(sc.embedded, features="c1_pcscore", reduction="fa", cols=viridis_pal()(100))
FeaturePlot(sc.embedded, features="c2_pcscore", reduction="fa", cols=viridis_pal()(100))
FeaturePlot(sc.embedded, features="c5_pcscore", reduction="fa", cols=viridis_pal()(100))
FeaturePlot(sc.embedded, features="c7_pcscore", reduction="fa", cols=viridis_pal()(100))
saveRDS(sc.embedded, "../data/seurat/genemodules.annotated.rds")

# directly compare module scores
FeatureScatter(sc, "c1_pcscore", "c2_pcscore", group.by="diffday")
FeatureScatter(sc, "c1_pcscore", "c5_pcscore", group.by="diffday")
FeatureScatter(sc, "c5_pcscore", "c9_pcscore", group.by="diffday")