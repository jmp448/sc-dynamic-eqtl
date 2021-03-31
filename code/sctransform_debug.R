library(Seurat)
library(tidyverse)
library(sctransform)

sc <- readRDS("../data/seurat.mitonormalized.rds")
sc <- RunPCA(sc)


sctransform::plot_model_pars(sc@assays$SCT@misc$vst.out, show_theta=T)

sc.sub <- sc[,sample(colnames(sc), 10000)]
sc.sub <- SCTransform(sc.sub)


sctransform::plot_model_pars(sc.sub@assays$SCT@misc$vst.out, show_theta=T)

# residual variance dist
ggplot(sc@assays$SCT@misc$vst.out$gene_attr, aes(residual_mean)) + 
  geom_histogram(binwidth = 0.1) +
  theme_classic()
ggplot(sc@assays$SCT@misc$vst.out$gene_attr, aes(residual_variance)) + 
  geom_histogram(binwidth = 0.1) + 
  geom_vline(xintercept = 1, color = "red") + 
  xlim(0, 10) +
  theme_classic()

# look at most variable genes
vst_out <- sc@assays$SCT@misc$vst.out
head(round(vst_out$gene_attr[order(-vst_out$gene_attr$residual_variance), ], 2), 
     22)