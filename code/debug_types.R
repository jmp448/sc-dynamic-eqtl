library(Seurat)
library(tidyverse)
library(Matrix)
library(Nebulosa)

sc <- readRDS("../data/seurat.annotated.rds")
DimPlot(sc)
plot_density(sc, "percent.mito")
PCAPlot(sc, group.by="type")

sc$pct.zero <- 1-colSums(sc@assays$SCT@counts!=0)/ncol(sc)
FeatureScatter(sc, "PC_2", "APOA2")
sc$low.mito <- sc$percent.mito<=2.5
DimPlot(sc, group.by="low.mito")
inc.cells <- read_tsv("../marcc_transfer2/graph_embedding.tsv") %>%
  .$X1
sc$inc <- colnames(sc) %in% inc.cells
DimPlot(sc, group.by="inc")

weirdos <- sc[,sc$type=="UNK"]

DimPlot(weirdos, group.by="diffday")
weirdos$pcdrive <- weirdos@reductions$pca@cell.embeddings[,2]>=50
PCAPlot(weirdos, group.by="pcdrive")
DimPlot(weirdos, group.by="pcdrive")
plot_density(weirdos, "APOA1")

super.weirdos <- weirdos[,weirdos$pcdrive==T]
table(super.weirdos$individual)

odd.sample <- sc[,sc$sample=="18855_11"]
table(odd.sample$type)
sc$odd <- sc$sample=="18489_11"
DimPlot(sc, group.by="odd")


weirdos2 <- weirdos %>%
  ScaleData %>%
  RunPCA %>%
  FindNeighbors %>%
  FindClusters %>%
  RunUMAP(dims=1:10)

DimPlot(weirdos2, group.by="diffday")
DimPlot(weirdos2, group.by="pcdrive")
FeaturePlot(weirdos2, features=c("nFeature_RNA", "nCount_RNA", "PRB.DBL"))
plot_density(weirdos2, "APOA1")
plot_density(weirdos2, "GATA4")
plot_density(weirdos2, "ISL1")
weirdos2$pcdrive <- colnames(weirdos2) %in% colnames(super.weirdos)

outliers_conservative <- as_tibble(weirdos2$seurat_clusters, rownames="cell") %>%
  filter(value %in% c(5)) %>%
  write_tsv("../data/outliers_conservative.tsv")
outliers_liberal <- as_tibble(weirdos2$seurat_clusters, rownames="cell") %>%
  filter(value %in% c(5,10,12)) %>%
  write_tsv("../data/outliers_liberal.tsv")

sc.cons <- readRDS("../data/seurat.conservative.rds")
sc2 <- readRDS("../data/seurat.liberal.rds")

sc.test <- FindVariableFeatures(sc)
"APOA2" %in% VariableFeatures(object=sc.test)

# look at cluster centroids
DimPlot(sc)
c0.umap.centroid <- sc@reductions$umap@cell.embeddings[sc$seurat_clusters==1,] %>% colMeans %>% as.matrix
c1.umap.centroid <- sc@reductions$umap@cell.embeddings[sc$seurat_clusters==2,] %>% colMeans %>% as.matrix
c2.umap.centroid <- sc@reductions$umap@cell.embeddings[sc$seurat_clusters==3,] %>% colMeans %>% as.matrix
c3.umap.centroid <- sc@reductions$umap@cell.embeddings[sc$seurat_clusters==4,] %>% colMeans %>% as.matrix
c4.umap.centroid <- sc@reductions$umap@cell.embeddings[sc$seurat_clusters==5,] %>% colMeans %>% as.matrix
c5.umap.centroid <- sc@reductions$umap@cell.embeddings[sc$seurat_clusters==6,] %>% colMeans %>% as.matrix
c6.umap.centroid <- sc@reductions$umap@cell.embeddings[sc$seurat_clusters==7,] %>% colMeans %>% as.matrix 
centroids <- t(cbind(c0.umap.centroid, c1.umap.centroid, c2.umap.centroid, c3.umap.centroid, c4.umap.centroid, c5.umap.centroid, c6.umap.centroid))
rownames(centroids) <- c("CF", "MESO", "IPSC", "PROG", "CMES", "CM", "UNK")
dist(centroids)

# now look at the subset data
sc.cons <- readRDS("../data/seurat.conservative.rds")
DimPlot(sc.cons)
npcs <- 100
c0.umap.centroid <- sc.cons@reductions$umap@cell.embeddings[sc.cons$seurat_clusters==0,1:npcs] %>% colMeans %>% as.matrix
c1.umap.centroid <- sc.cons@reductions$umap@cell.embeddings[sc.cons$seurat_clusters==1,1:npcs] %>% colMeans %>% as.matrix
c2.umap.centroid <- sc.cons@reductions$umap@cell.embeddings[sc.cons$seurat_clusters==2,1:npcs] %>% colMeans %>% as.matrix
c3.umap.centroid <- sc.cons@reductions$umap@cell.embeddings[sc.cons$seurat_clusters==3,1:npcs] %>% colMeans %>% as.matrix
c4.umap.centroid <- sc.cons@reductions$umap@cell.embeddings[sc.cons$seurat_clusters==4,1:npcs] %>% colMeans %>% as.matrix
c5.umap.centroid <- sc.cons@reductions$umap@cell.embeddings[sc.cons$seurat_clusters==5,1:npcs] %>% colMeans %>% as.matrix
c6.umap.centroid <- sc.cons@reductions$umap@cell.embeddings[sc.cons$seurat_clusters==6,1:npcs] %>% colMeans %>% as.matrix 
c7.umap.centroid <- sc.cons@reductions$umap@cell.embeddings[sc.cons$seurat_clusters==7,1:npcs] %>% colMeans %>% as.matrix 
c8.umap.centroid <- sc.cons@reductions$umap@cell.embeddings[sc.cons$seurat_clusters==8,1:npcs] %>% colMeans %>% as.matrix 
c9.umap.centroid <- sc.cons@reductions$umap@cell.embeddings[sc.cons$seurat_clusters==9,1:npcs] %>% colMeans %>% as.matrix 
centroids <- t(cbind(c0.umap.centroid, c1.umap.centroid, c2.umap.centroid, c3.umap.centroid, c4.umap.centroid, 
                     c5.umap.centroid, c6.umap.centroid, c7.umap.centroid, c8.umap.centroid, c9.umap.centroid))
rownames(centroids) <- c("MESO", "PROG/CF", "IPSC", "CM", "CMES", "UNK", "BRIDGE", "CF2", "INVIS1", "INVIS2")
dist(centroids)
plot(hclust(dist(centroids)))

c0.umap.centroid <- sc.cons@reductions$umap@cell.embeddings[sc.cons$seurat_clusters==0,] %>% colMeans %>% as.matrix
c1.umap.centroid <- sc.cons@reductions$umap@cell.embeddings[sc.cons$seurat_clusters==1,] %>% colMeans %>% as.matrix
c2.umap.centroid <- sc.cons@reductions$umap@cell.embeddings[sc.cons$seurat_clusters==2,] %>% colMeans %>% as.matrix
c3.umap.centroid <- sc.cons@reductions$umap@cell.embeddings[sc.cons$seurat_clusters==3,] %>% colMeans %>% as.matrix
c4.umap.centroid <- sc.cons@reductions$umap@cell.embeddings[sc.cons$seurat_clusters==4,] %>% colMeans %>% as.matrix
c5.umap.centroid <- sc.cons@reductions$umap@cell.embeddings[sc.cons$seurat_clusters==5,] %>% colMeans %>% as.matrix
c6.umap.centroid <- sc.cons@reductions$umap@cell.embeddings[sc.cons$seurat_clusters==6,] %>% colMeans %>% as.matrix 
c7.umap.centroid <- sc.cons@reductions$umap@cell.embeddings[sc.cons$seurat_clusters==7,] %>% colMeans %>% as.matrix 
c8.umap.centroid <- sc.cons@reductions$umap@cell.embeddings[sc.cons$seurat_clusters==8,] %>% colMeans %>% as.matrix 
c9.umap.centroid <- sc.cons@reductions$umap@cell.embeddings[sc.cons$seurat_clusters==9,] %>% colMeans %>% as.matrix 
centroids.umap <- t(cbind(c0.umap.centroid, c1.umap.centroid, c2.umap.centroid, c3.umap.centroid, c4.umap.centroid, 
                     c5.umap.centroid, c6.umap.centroid, c7.umap.centroid, c8.umap.centroid, c9.umap.centroid))
rownames(centroids.umap) <- c("MESO", "PROG/CF", "IPSC", "CM", "CMES", "UNK", "BRIDGE", "CF2", "INVIS1", "INVIS2")
dist(centroids.umap)
plot(hclust(dist(centroids.umap)))

# can we get the force atlas embedding into seurat?
fa <- read_tsv("../marcc_transfer2/graph_embedding.tsv") %>%
  rename(FA_1=`0`, FA_2=`1`) %>% column_to_rownames("X1") %>%
  as.matrix
test <- sc[,colnames(sc) %in% rownames(fa)]
fa <- fa[colnames(test),]
test[["fa"]] <- CreateDimReducObject(embeddings = fa, key = "FA_", assay = DefaultAssay(test))
DimPlot(test, reduction="fa", pt.size=0.001)
  

sc.cons <- readRDS("../data/seurat.conservative.rds")
PCAPlot(sc.cons, c(1,4))
PCASigGenes(sc.cons, pcs.use=1:2)

ProjectDim(sc.cons, reduction="pca", dims.print=1:10)
FeatureScatter(sc.cons, "PC_5", "S.Score")

# compare the cells whose assignment changes to the original members of that group
leiden <- read_tsv("../marcc_transfer2/leiden.tsv")
leiden_md <- leiden$leiden
names(leiden_md) <- leiden$X1
sc <- AddMetaData(sc, leiden_md, col.name="leiden")
DimPlot(sc, group.by="leiden")

mids <- sc[,sc$leiden %in% c("3", "0")]
DimPlot(mids, group.by="leiden")
mids <- mids %>%
  ScaleData %>%
  RunPCA %>%
mids <- mids %>%
  FindNeighbors %>%
  FindClusters %>%
  RunUMAP(dims=1:10)
DimPlot(mids, group.by="type")
FindMarkers(mids, ident.1=4)

log.mean <- sc@assays$SCT@data %>% apply(1, mean)
log.var <- sc@assays$SCT@data %>% apply(1, var)
log.info <- tibble("gene"=rownames(sc@assays$SCT@data), "mean"=log.mean, "var"=log.var)

res.mean <- sc@assays$SCT@scale.data %>% apply(1, mean)
res.var <- sc@assays$SCT@scale.data %>% apply(1, var)
res.info <- tibble("gene"=rownames(sc@assays$SCT@scale.data), "mean"=res.mean, "var"=res.var)

log.info <- 
ggplot(log.info, aes(x=mean)) + geom_histogram()
ggplot(log.info, aes(x=var)) + geom_histogram(bins=100)
ggplot(res.info, aes(x=mean)) + geom_histogram()
ggplot(res.info, aes(x=var)) + geom_histogram(bins=100) +
  geom_vline(xintercept=1, color="red") +
  xlim(0, 10) +
  theme_classic()
test <- filter(res.info, substr(res.info$gene, 1, 3)=="MT-")
head(test)

gene.mean <- sc@assays$SCT@meta.features$sct.gmean
res.var <- sc@assays$SCT@meta.features$sct.residual_variance
res.mean <- sc@assays$SCT@meta.features$sct.residual_mean
sctrans <- as_tibble(sc@assays$SCT@meta.features)
sctrans$gene <- rownames(sc)
ggplot(sctrans, aes(x=gene.mean, y=res.var)) +
  geom_point() +
  scale_x_continuous(trans='log10')

apoa1 <- res.info %>% filter(gene %in% c("APOA1", "APOA2", "AFP"))
marks <- res.info %>% filter(gene %in% c("TNNT2", "COL3A1", "L1TD1", "MIXL1", "MESP1"))