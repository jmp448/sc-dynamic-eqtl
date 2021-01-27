library(tidyverse)
library(edgeR)
library(Seurat)
library(Matrix)
library(Matrix.utils)
library(lars)
library(Nebulosa)
library(glmnet)
set.seed(2021)
source("helpers.R")
sc <- readRDS("/project2/gilad/jpopp/sc3/data/mitofilt_singlets.full-clean.annotated.cc.rds")

# for now, using day 11 cell type proportions
ct.proportions <- read_tsv("../data/cibersort/outputs/pseudobulk_labels.true.tsv")
ct.proportions <- ct.proportions %>% filter(grepl("_11", sample)) 
ct.proportions$ind <- paste0("NA", str_sub(ct.proportions$sample, 1, 5))
props <- sc$individual %>% as_tibble(rownames="cell") %>% rename(ind=value) %>% left_join(ct.proportions, by="ind")
sc$cm.day11 <- props$CM
sc$epdc.day11 <- props$EPDC
sc$emt.day11 <- props$EMT

# subset to day 0 cells
sc.ipsc <- sc[,(sc$diffday=="Day 0" & sc$type=="iPSC")] 
sc.ipsc <- sc.ipsc %>% FindVariableFeatures %>% 
  ScaleData(vars.to.regress=c("G2M.Score", "S.Score")) %>% 
  RunPCA %>% FindNeighbors %>% FindClusters(res=0.5) %>% RunUMAP(dims=1:10)
saveRDS(sc.ipsc, "../data/seurat/ipsc.rds")

DimPlot(sc.ipsc, group.by="individual")
plot_density(sc.ipsc, "UTF1")
plot_density(sc.ipsc, "TAC3")
plot_density(sc.ipsc, "cm.day11")
plot_density(sc.ipsc, "epdc.day11")
plot_density(sc.ipsc, "emt.day11")

# correlation between pseudobulk and CM fraction
pseudobulk <- read_tsv("../data/pseudobulk_cpm.day.overlap.tsv")
pseudobulk_labels <- read_tsv("../data/cibersort/outputs/pseudobulk_labels.true.tsv")

cm_labels <- pseudobulk_labels %>% filter(grepl("_11", sample)) %>% 
  mutate(sample=str_sub(sample, 1, 5)) %>% select(c(sample, CM))
ipsc_exp <- pseudobulk %>% filter(gene %in% c("UTF1", "TAC3")) %>%
  select(c(gene, contains("_0"))) %>% column_to_rownames("gene") %>%
  t %>% as_tibble(rownames="sample") %>% mutate(sample=str_sub(sample, 1, 5))
predictors <- inner_join(cm_labels, ipsc_exp, by="sample")
ggplot(predictors, aes(x=UTF1, y=CM)) + geom_point()
ggplot(predictors, aes(x=TAC3, y=CM)) + geom_point()

# scPCs and CM fraction
FeaturePlot(sc.ipsc, "cm.day11", )
samples.ranked <- cm_labels %>% arrange(CM) %>% .$sample
sc.ipsc$individual <- sc.ipsc$individual %>% str_sub(3) %>% factor(levels=samples.ranked)
VlnPlot(sc.ipsc, "PC_1", pt.size=0, group.by="individual")
VlnPlot(sc.ipsc, "PC_2", pt.size=0, group.by="individual")
VlnPlot(sc.ipsc, "PC_3", pt.size=0, group.by="individual")
VlnPlot(sc.ipsc, "PC_4", pt.size=0, group.by="individual")
VlnPlot(sc.ipsc, "PC_5", pt.size=0, group.by="individual")
cm_labels <- cm_labels %>% arrange(CM) %>% mutate(sample=factor(sample, levels=samples.ranked))
ggplot(cm_labels, aes(x=sample, y=CM, fill=sample)) + 
  geom_bar(stat="identity") +
  theme(axis.text.x=element_text(angle=45))


# subset to day 1 meso cells
sc$d1 <- (sc$diffday == "Day 1" & sc$type == "meso")
DimPlot(sc, group.by="d1")
sc.meso <- sc[,(sc$diffday=="Day 1" & sc$type=="meso")] 
sc.meso <- sc.meso %>% FindVariableFeatures %>% 
  ScaleData(vars.to.regress=c("G2M.Score", "S.Score")) %>% 
  RunPCA %>% FindNeighbors %>% FindClusters(res=0.5) %>% RunUMAP(dims=1:10)
saveRDS(sc.meso, "../data/seurat/meso.rds")

plot_density(sc.meso, "cm.day11")
plot_density(sc.meso, "epdc.day11")
plot_density(sc.meso, "emt.day11")
