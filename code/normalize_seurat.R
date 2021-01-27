library(Seurat)
library(tidyverse)

sc <- readRDS("../data/seurat.filtered.rds")

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
sc <- CellCycleScoring(sc, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
sc$CC.Difference <- sc$S.Score - sc$G2M.Score

options(future.globals.maxSize=1*1000^3)
sc <- SCTransform(sc, vars.to.regress=c("CC.Difference"))
saveRDS(sc, "../data/seurat.normalized.rds")