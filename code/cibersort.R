library(tidyverse)
library(edgeR)
library(Seurat)
library(Matrix)
library(Matrix.utils)
library(corrplot)
set.seed(2021)

# load bulk data
bulk <- read_tsv("../data/cibersort/inputs/bulk_tpm.full.tsv")

# load single cell data
sc <- readRDS("/project2/gilad/jpopp/sc3/data/mitofilt_singlets.full-clean.annotated.cc.rds")

# train/test split (60/40) for assessment of cell type deconvolution
cell.ids <- as_tibble(sc$type, rownames="cell")
test.cells <- cell.ids %>% group_by(value) %>% sample_frac(0.4) %>% .$cell
train.cells <- setdiff(cell.ids$cell, test.cells)
sc.test <- sc[,test.cells]
sc.train <- sc[,train.cells]

# get differentially expressed genes to use for deconvolution
Idents(sc.train) <- sc.train$type
ct.markers <- FindAllMarkers(sc.train, only.pos=T)
sig.genes <- ct.markers$gene

# subset sig.genes to genes that are also measured in bulk
bulk.genes <- bulk$gene
sig.genes <- intersect(sig.genes, bulk.genes)

write_tsv(tibble(sig.genes), "../data/cibersort/sig.genes.tsv")

# SIGNATURE MATRIX
# aggregate training data into pseudobulk by type
sig.matrix <- sc.train@assays[["SCT"]]@counts %>% t %>% 
  aggregate.Matrix(sc.train$type, fun="sum") %>% 
  t %>% as_tibble(rownames="gene")
# normalize to get CPM data
sig.matrix <- sig.matrix %>% select(!gene) %>% DGEList %>% calcNormFactors(method="TMMwsp") %>% 
  cpm %>% as_tibble %>% mutate(gene=sig.matrix$gene) %>%
  relocate(gene)
# filter to marker genes
sig.matrix <- filter(sig.matrix, gene %in% sig.genes)

write_tsv(sig.matrix, "../data/cibersort/inputs/signature_matrix.tsv")

# ASSESSMENT 
# aggregate by (individual/collection time point) - separates replicates
sc.test$sample <- paste( sc.test$individual,sc.test$orig.ident,sep="_") %>% 
  str_replace("NA", "")
pb.test <- sc.test@assays[["SCT"]]@counts %>% t %>% 
  aggregate.Matrix(sc.test$sample, fun="sum") %>% 
  t %>% as_tibble(rownames="gene")
# normalize to get CPM data
pb.test <- pb.test %>% select(!gene) %>% DGEList %>% calcNormFactors(method="TMMwsp") %>% 
  cpm %>% as_tibble %>% mutate(gene=pb.test$gene) %>%
  relocate(gene)
# subset to genes evaluated for deconvolution
pb.test <- pb.test %>% filter(gene %in% sig.genes)

write_tsv(pb.test, "../data/cibersort/inputs/pseudobulk_cpm.indcol.test.sig.tsv")

# ground truth 
true.props <- tibble("sample"=sc.test$sample, "type"=sc.test$type) %>%
  mutate(n=1) %>%
  group_by(sample, type) %>%
  summarize(count=sum(n)) %>% ungroup
sampsums <- true.props %>% group_by(sample) %>%
  summarize(sampsum=sum(count))
true.props <- left_join(true.props, sampsums, by="sample") %>%
  mutate(frac=count/sampsum) %>%
  select(!c(count, sampsum)) %>%
  spread(key=type, value=frac) %>%
  mutate_all(~replace_na(.,0))

write_tsv(true.props, "../results/cibersort/pseudobulk_fracs.indcol.test.true.tsv")
  
# aggregate by (individual/development time point) - merges replicates
sc.test$sample <- paste(sc.test$individual, sc.test$diffday, sep="_") %>%
  str_replace("NA", "") %>% str_replace("Day ", "")
pb.test <- sc.test@assays[["SCT"]]@counts %>% t %>% 
  aggregate.Matrix(sc.test$sample, fun="sum") %>% 
  t %>% as_tibble(rownames="gene")
# normalize to get CPM data
pb.test <- pb.test %>% select(!gene) %>% DGEList %>% calcNormFactors(method="TMMwsp") %>% 
  cpm %>% as_tibble %>% mutate(gene=pb.test$gene) %>%
  relocate(gene)
# subset to genes evaluated for deconvolution
pb.test <- pb.test %>% filter(gene %in% sig.genes)

write_tsv(pb.test, "../data/cibersort/inputs/pseudobulk_cpm.indday.test.sig.tsv")

# ground truth 
true.props <- tibble("sample"=sc.test$sample, "type"=sc.test$type) %>%
  mutate(n=1) %>%
  group_by(sample, type) %>%
  summarize(count=sum(n)) %>% ungroup
sampsums <- true.props %>% group_by(sample) %>%
  summarize(sampsum=sum(count))
true.props <- left_join(true.props, sampsums, by="sample") %>%
  mutate(frac=count/sampsum) %>%
  select(!c(count, sampsum)) %>%
  spread(key=type, value=frac) %>%
  mutate_all(~replace_na(.,0))

write_tsv(true.props, "../results/cibersort/pseudobulk_fracs.indday.test.true.tsv")

# ground truth for full pseudobulk data
sc$sample <- paste(sc$individual, sc$diffday, sep="_") %>%
  str_replace("NA", "") %>% str_replace("Day ", "")
true.props <- tibble("sample"=sc$sample, "type"=sc$type) %>%
  mutate(n=1) %>%
  group_by(sample, type) %>%
  summarize(count=sum(n)) %>% ungroup
sampsums <- true.props %>% group_by(sample) %>%
  summarize(sampsum=sum(count))
true.props <- left_join(true.props, sampsums, by="sample") %>%
  mutate(frac=count/sampsum) %>%
  select(!c(count, sampsum)) %>%
  spread(key=type, value=frac) %>%
  mutate_all(~replace_na(.,0))

write_tsv(true.props, "../results/cibersort/pseudobulk_fracs.indday.full.true.tsv")

sc$sample <- paste(sc$individual,sc$orig.ident,sep="_") %>% 
  str_replace("NA", "")
true.props <- tibble("sample"=sc$sample, "type"=sc$type) %>%
  mutate(n=1) %>%
  group_by(sample, type) %>%
  summarize(count=sum(n)) %>% ungroup
sampsums <- true.props %>% group_by(sample) %>%
  summarize(sampsum=sum(count))
true.props <- left_join(true.props, sampsums, by="sample") %>%
  mutate(frac=count/sampsum) %>%
  select(!c(count, sampsum)) %>%
  spread(key=type, value=frac) %>%
  mutate_all(~replace_na(.,0))

write_tsv(true.props, "../results/cibersort/pseudobulk_fracs.indcol.full.true.tsv")

# subset bulk data for deconvolution
bulk <- read_tsv("../data/bulk_tpm.full.tsv") %>%
  filter(gene %in% sig.genes$sig.genes)
write_tsv(bulk, "../data/cibersort/inputs/bulk_tpm.sig.tsv")