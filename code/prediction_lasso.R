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

# load pseudobulk data and cell type proportions
pseudobulk <- read_tsv("../data/pseudobulk_counts.day.overlap.tsv")
pseudobulk_labels <- read_tsv("../data/cibersort/outputs/pseudobulk_labels.inferred.tsv")

# match pseudobulk column order to label row order
pseudobulk <- pseudobulk %>% relocate(c("gene", pseudobulk_labels$sample))

# LASSO regression of cell type proportions onto expression profiles
d.design <- "11"
d.response <- "11"
pb <- pseudobulk[,str_extract(colnames(pseudobulk), "[^_]+$")==d.design] %>% DGEList %>% 
  calcNormFactors(method="TMMwsp") %>% cpm(log=T) %>%
  apply(1, center.scale) %>% t %>% as_tibble %>% mutate(gene=pseudobulk$gene) %>%
  relocate(gene) %>% arrange(gene)
pb_labels <- pseudobulk_labels %>% filter(grepl(paste0("_", d.response), sample))

design <- pb %>% select(!gene) %>% as.matrix %>% t %>% `colnames<-`(pseudobulk$gene)
response <- pb_labels %>% .$CM %>% `names<-`(pb_labels$sample)
rownames(design) <- str_extract(rownames(design), "[^_]+")
names(response) <- str_extract(names(response), "[^_]+")
response <- response[rownames(design)]
assert_that(sum(names(response)!=rownames(design))==0)

lasso.all <- cv.glmnet(design, response, family="gaussian", alpha=1, nfolds=19)
lambduh <- lasso.all$lambda.min
cvm <- lasso.all$cvm[lasso.all$lambda==lambduh]
lasso.best <- glmnet(design, response, alpha=1, family="gaussian", lambda=lambduh)
coefs <- as_tibble(lasso.best$beta, rownames="gene") %>% rename(beta=s0)
loaded.coefs <- coefs %>% filter(beta!=0) %>% arrange(gene)

# get predictors, and predictor genes (to measure R2) - pseudobulk
B.pb <- loaded.coefs %>% select(!gene) %>% as.matrix %>% `rownames<-`(loaded.coefs$gene)
X.pb <- pb %>% filter(gene %in% loaded.coefs$gene) %>% arrange(gene) %>%
  select(!gene) %>% as.matrix %>% t %>% `colnames<-`(loaded.coefs$gene)
CM.hat <- X.pb %*% B.pb
comp <- CM.hat %>% as_tibble(rownames="sample") %>% rename(CM.hat=beta) %>% inner_join(select(pb_labels, c(sample, CM)), by="sample")
ggplot(comp, aes(x=CM, y=CM.hat)) + geom_point()

# get predictions for bulk
B.bulk <- B.pb
X.bulk <- 

# look at the expression patterns of these predictive genes
plot_density(sc, "MCFD2")
pb %>% filter(gene=="MCFD2") %>% select(!gene) %>% as.numeric %>% 
  as_tibble %>% rename(WNK2=value) %>% mutate(sample=colnames(pb)[-c(1)]) %>%
  inner_join(pb_labels, by="sample") %>% select(WNK2, CM) %>%
  ggplot(aes(x=WNK2, y=CM)) + geom_point()

