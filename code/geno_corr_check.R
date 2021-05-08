library(tidyverse)
library(limma)
library(vroom)
library(qvalue)
library(Seurat)
library(gplots)

dataset <- "bulk"
agg <- "day"
all.tests <- vroom(paste0("../data/", dataset, "/", agg, "/filtered_tests.50k.tsv"))
dyn.qtls <- read_tsv(paste0("../results/eqtl_dynamic/linear_dQTL/", dataset, "/", agg, "/50k-5clpcs-0pcs-notypes.tophits.tsv")) %>%
  filter(qval.unadj<=0.05) %>%
  mutate(gv=paste(gene, snp, sep="--"))

mafs <- vroom("../data/geno_maf.tsv")
dqtl.top200 <- dyn.qtls %>% 
  arrange(qval.unadj) %>% 
  slice_head(n=200) %>% 
  select(gene, snp, gv) %>%
  mutate(group="dynqtl") %>%
  mutate("ct_regressed"=F) %>%
  left_join(mafs, by="snp") %>%
  mutate(maf.bin=floor(maf/0.05))
bin.counts <- count(dqtl.top200, maf.bin)
dqtl.top200 <- select(dqtl.top200, snp) 
# match MAF for some background variants
matched.top200 <- all.tests %>%
  left_join(mafs, by="snp")%>%
  mutate(maf.bin=floor(maf/0.05)) %>%
  group_by(maf.bin) %>%
  nest() %>%
  ungroup() %>%
  left_join(bin.counts, by="maf.bin") %>%
  replace_na(list(n=0)) %>%
  mutate(samp=map2(data, n, sample_n)) %>%
  select(!data) %>%
  unnest(samp) %>%
  select(snp)

# get a (nsnps x ncellline) matrix of each individual's genotype at each snp
genotypes <- vroom("../data/genotypes.tsv")
dqtl.cor <- genotypes %>% 
  filter(snp %in% dqtl.top200$snp) %>%
  column_to_rownames("snp") %>%
  cor
matched.cor <- genotypes %>%
  filter(snp %in% matched.top200$snp) %>%
  column_to_rownames("snp") %>%
  cor 

corr.cols <- colorRampPalette(brewer.pal(11, "RdBu"))(25)

png(paste0('../figs/supp/dqtl_corr_', dataset, '_', agg, '_hits.png'), width=800, height=600)
heatmap.2(dqtl.cor, Rowv=T, Colv=T, dendrogram='none', trace='none', 
          col=rev(corr.cols), breaks=seq(-1, 1, length.out=26),
          main="Bulk, 16 Days, Dynamic eQTL Variants")
dev.off()

png(paste0('../figs/supp/dqtl_corr_', dataset, '_', agg, '_matched.png'), width=800, height=600)
heatmap.2(matched.cor, Rowv=T, Colv=T, dendrogram='none', trace='none', 
          col=rev(corr.cols), breaks=seq(-1, 1, length.out=26),
          main="Bulk, 16 Days, Background Variants")
dev.off()

### interaction eqtl version
dataset <- "bulk"
ct <- "CF"
notypes.qtls <- read_tsv(paste0("../results/eqtl_dynamic/ieQTL/bulk/CF/50k-5clpcs-0pcs-notypes.tophits.tsv")) %>%
  filter(qval.unadj<=0.05) %>%
  mutate(gv=paste(gene, snp, sep="--"))
regtypes.qtls <- read_tsv(paste0("../results/eqtl_dynamic/ieQTL/bulk/CF/50k-5clpcs-0pcs-regtypes.tophits.tsv")) %>%
  filter(qval.unadj<=0.05) %>%
  filter(!snp %in% notypes.qtls$snp) %>%
  mutate(gv=paste(gene, snp, sep="--"))
pc5.qtls <- read_tsv(paste0("../results/eqtl_dynamic/ieQTL/bulk/CF/50k-5clpcs-5pcs-notypes.tophits.tsv")) %>%
  filter(qval.unadj<=0.05) %>%
  filter(!snp %in% notypes.qtls$snp) %>%
  mutate(gv=paste(gene, snp, sep="--"))
pc10.qtls <- read_tsv(paste0("../results/eqtl_dynamic/ieQTL/bulk/CF/50k-5clpcs-10pcs-notypes.tophits.tsv")) %>%
  filter(qval.unadj<=0.05) %>%
  filter(!snp %in% pc5.qtls$snp) %>%
  mutate(gv=paste(gene, snp, sep="--"))
pc20.qtls <- read_tsv(paste0("../results/eqtl_dynamic/ieQTL/bulk/CF/50k-5clpcs-20pcs-notypes.tophits.tsv")) %>%
  filter(qval.unadj<=0.05) %>%
  filter(!snp %in% pc10.qtls$snp) %>%
  mutate(gv=paste(gene, snp, sep="--"))
pc30.qtls <- read_tsv(paste0("../results/eqtl_dynamic/ieQTL/bulk/CF/50k-5clpcs-30pcs-notypes.tophits.tsv")) %>%
  filter(qval.unadj<=0.05) %>%
  filter(!snp %in% pc20.qtls$snp) %>%
  mutate(gv=paste(gene, snp, sep="--"))

notypes.top200 <- notypes.qtls %>% 
  arrange(qval.unadj) %>% 
  slice_head(n=200) %>% 
  select(snp) 

regtypes.top200 <- regtypes.qtls %>% 
  arrange(qval.unadj) %>% 
  slice_head(n=200) %>% 
  select(snp) 

pc5.top200 <- pc5.qtls %>% 
  arrange(qval.unadj) %>% 
  slice_head(n=200) %>% 
  select(snp) 

pc10.top200 <- pc10.qtls %>% 
  arrange(qval.unadj) %>% 
  slice_head(n=200) %>% 
  select(snp) 

pc20.top200 <- pc20.qtls %>% 
  arrange(qval.unadj) %>% 
  slice_head(n=200) %>% 
  select(snp) 

pc30.top200 <- pc30.qtls %>% 
  arrange(qval.unadj) %>% 
  slice_head(n=200) %>% 
  select(snp) 


# get a (nsnps x ncellline) matrix of each individual's genotype at each snp
notypes.cor <- genotypes %>% 
  filter(snp %in% notypes.top200$snp) %>%
  column_to_rownames("snp") %>%
  cor
regtypes.cor <- genotypes %>% 
  filter(snp %in% regtypes.top200$snp) %>%
  column_to_rownames("snp") %>%
  cor
pc5.cor <- genotypes %>% 
  filter(snp %in% pc5.top200$snp) %>%
  column_to_rownames("snp") %>%
  cor
pc10.cor <- genotypes %>% 
  filter(snp %in% pc10.top200$snp) %>%
  column_to_rownames("snp") %>%
  cor
pc20.cor <- genotypes %>% 
  filter(snp %in% pc20.top200$snp) %>%
  column_to_rownames("snp") %>%
  cor
pc30.cor <- genotypes %>% 
  filter(snp %in% pc30.top200$snp) %>%
  column_to_rownames("snp") %>%
  cor

corr.cols <- colorRampPalette(brewer.pal(11, "RdBu"))(25)

png(paste0('../figs/supp/ieqtl_corr_notypes.png'), width=800, height=600)
heatmap.2(notypes.cor, Rowv=T, Colv=T, dendrogram='none', trace='none', 
          col=rev(corr.cols), breaks=seq(-1, 1, length.out=26),
          main="No Additional Covariates")
dev.off()

png(paste0('../figs/supp/ieqtl_corr_regtypes.png'), width=800, height=600)
heatmap.2(regtypes.cor, Rowv=T, Colv=T, dendrogram='none', trace='none', 
          col=rev(corr.cols), breaks=seq(-1, 1, length.out=26),
          main="All Cell Type Covariates")
dev.off()

png(paste0('../figs/supp/ieqtl_corr_5pc.png'), width=800, height=600)
heatmap.2(pc5.cor, Rowv=T, Colv=T, dendrogram='none', trace='none', 
          col=rev(corr.cols), breaks=seq(-1, 1, length.out=26),
          main="5 Sample PCs")
dev.off()

png(paste0('../figs/supp/ieqtl_corr_10pc.png'), width=800, height=600)
heatmap.2(pc10.cor, Rowv=T, Colv=T, dendrogram='none', trace='none', 
          col=rev(corr.cols), breaks=seq(-1, 1, length.out=26),
          main="10 Sample PCs")
dev.off()

png(paste0('../figs/supp/ieqtl_corr_20pc.png'), width=800, height=600)
heatmap.2(pc20.cor, Rowv=T, Colv=T, dendrogram='none', trace='none', 
          col=rev(corr.cols), breaks=seq(-1, 1, length.out=26),
          main="20 Sample PCs")
dev.off()

png(paste0('../figs/supp/ieqtl_corr_30pc.png'), width=800, height=600)
heatmap.2(pc30.cor, Rowv=T, Colv=T, dendrogram='none', trace='none', 
          col=rev(corr.cols), breaks=seq(-1, 1, length.out=26),
          main="30 Sample PCs")
dev.off()
