library(assertthat)
library(tidyverse)
library(edgeR)
library(Seurat)
library(Matrix)
library(Matrix.utils)
library(corrplot)
set.seed(2021)
source("/project2/gilad/jpopp/sc-dynamic-eqtl/code/cell_line_pca.R")

# BULK 
# load the bulk data, convert to RPKM (to compare to Ben's) and TPM
bulk <- read_delim("/project2/gilad/bstrober/ipsc_differentiation_19_lines/preprocess/processed_total_expression/raw_counts.txt",
                   delim=" ")
gene_len <- read_tsv("/project2/gilad/bstrober/ipsc_differentiation_19_lines/preprocess/processed_total_expression/gene_lengths.txt") %>%
  rename(`Gene_id`=GeneID)
full_data <- bulk %>% inner_join(gene_len, by="Gene_id")
counts <- full_data %>% select(!c(`Gene_id`, Length))
lengths <- full_data[,"Length"]
rpk <- counts / lengths$Length
bulk.cpm <- DGEList(counts=counts) %>% calcNormFactors %>% cpm %>%
  as_tibble %>% mutate(gene=full_data$Gene_id) %>% relocate(gene)
bulk.tpm <- DGEList(counts=rpk) %>% calcNormFactors %>% cpm %>%
  as_tibble %>% mutate(gene=full_data$Gene_id) %>% relocate(gene)
bulk.rpkm <- DGEList(counts=counts, genes=lengths) %>% calcNormFactors %>% rpkm %>%
  as_tibble %>% mutate(gene=full_data$Gene_id) %>% relocate(gene)
bulk.logtpm <- DGEList(counts=rpk) %>% calcNormFactors %>% cpm(log=T) %>%
  as_tibble %>% mutate(gene=full_data$Gene_id) %>% relocate(gene)

# swap ENSG gene names for HGNC symbols
translator <- read_tsv("/project2/gilad/jpopp/cm_eqtl/data/gencode.v26.annotation.gene.txt") %>%
  mutate(gene=str_extract(gene_id, "[^.]+")) %>%
  filter(!duplicated(gene)) %>% select(gene, gene_name)
bulk.raw <- full_data %>% select(!Length) %>% rename(gene=`Gene_id`) %>%
  inner_join(translator, by="gene") %>% select(!gene) %>%
  rename(gene=gene_name) %>% relocate(gene) %>% arrange(gene)
bulk.cpm <- bulk.cpm %>% inner_join(translator, by="gene") %>% select(!gene) %>% 
  rename(gene=gene_name) %>% relocate(gene) %>% arrange(gene)
bulk.rpkm <- bulk.rpkm %>% inner_join(translator, by="gene") %>% select(!gene) %>% 
  rename(gene=gene_name) %>% relocate(gene) %>% arrange(gene)
bulk.tpm <- bulk.tpm %>% inner_join(translator, by="gene") %>% select(!gene) %>% 
  rename(gene=gene_name) %>% relocate(gene) %>% arrange(gene)
bulk.logtpm <- bulk.logtpm %>% inner_join(translator, by="gene") %>% select(!gene) %>% 
  rename(gene=gene_name) %>% relocate(gene) %>% arrange(gene)

assert_that(sum(colnames(bulk.raw) != colnames(bulk.cpm))==0)
assert_that(sum(colnames(bulk.raw) != colnames(bulk.rpkm))==0)
assert_that(sum(colnames(bulk.raw) != colnames(bulk.tpm))==0)
assert_that(sum(colnames(bulk.raw) != colnames(bulk.logtpm))==0)

write_tsv(bulk.raw, "../data/bulk_counts.full.tsv")
write_tsv(bulk.cpm, "../data/bulk_cpm.full.tsv")
write_tsv(bulk.rpkm, "../data/bulk_rpkm.full.tsv")
write_tsv(bulk.tpm, "../data/bulk_tpm.full.tsv")
write_tsv(bulk.logtpm, "../data/bulk_logtpm.full.tsv")

# PSEUDOBULK
sc <- readRDS("/project2/gilad/jpopp/sc3/data/mitofilt_singlets.full-clean.annotated.cc.rds")
counts <- sc@assays[["SCT"]]@counts %>% t

sc$type.ind <- paste(sc$individual, sc$type, sep="_") %>% str_replace("NA", "")
counts.type <- counts %>% aggregate.Matrix(sc$type.ind, fun="sum") %>% t %>% 
  as_tibble(rownames="gene")

sc$day.ind <- paste(sc$individual, sc$diffday, sep="_") %>% str_replace("NA", "") %>%
  str_replace("Day ", "")
counts.day <- counts %>% aggregate.Matrix(sc$day.ind, fun="sum") %>% t %>%
  as_tibble(rownames="gene")
day.drop <- counts.day %>% select(!gene) %>% colSums %>% 
  as_tibble(rownames="sample") %>% filter(value <= 5000) %>% .$sample
counts.day <- counts.day %>% select(-c(day.drop))

sc$col.ind <- paste(sc$individual, sc$orig.ident, sep="_") %>% 
  str_replace("NA", "")
counts.col <- counts %>% aggregate.Matrix(sc$col.ind, fun="sum") %>% t %>% 
  as_tibble(rownames="gene")
col.drop <- counts.col %>% select(!gene) %>% colSums %>% 
  as_tibble(rownames="sample") %>% filter(value <= 5000) %>% .$sample
counts.col <- counts.col %>% select(-c(col.drop))

# bin by pseudotime for cardiomyocyte and EPDC
cm.bins <- as_tibble(sc$CM_PC1, rownames="cell") %>% 
  drop_na %>% rename(t=value) %>% arrange(t) %>% 
  rowid_to_column("i") %>% mutate(bin=1+floor(15*i/length(i)), .keep="unused") 
cm.bins$bin[nrow(cm.bins)] = 15 # so the last cell doesn't get its own bin
cm.bins <- sc$individual %>% as_tibble(rownames="cell") %>%
  mutate(ind=str_replace(value, "NA", ""), .keep="unused") %>%
  right_join(cm.bins, by="cell") %>%
  mutate(bin.ind=paste(ind, bin, sep="_"))
cm.medians <- cm.bins %>%
  group_by(bin.ind) %>%
  summarize(t=median(t))
write_tsv(cm.medians, "../results/eqtl_dynamic/linear_dQTL/pseudobulk/cmbin/15bin_medians.tsv")

epdc.bins <- as_tibble(sc$EPDC_PC1, rownames="cell") %>% 
  drop_na %>% rename(t=value) %>% arrange(t) %>% 
  rowid_to_column("i") %>% mutate(bin=1+floor(15*i/length(i)), .keep="unused") 
epdc.bins$bin[nrow(epdc.bins)] = 15 # so the last cell doesn't get its own bin
epdc.bins <- sc$individual %>% as_tibble(rownames="cell") %>%
  mutate(ind=str_replace(value, "NA", ""), .keep="unused") %>%
  right_join(epdc.bins, by="cell") %>%
  mutate(bin.ind=paste(ind, bin, sep="_"))
epdc.medians <- epdc.bins %>%
  group_by(bin.ind) %>%
  summarize(t=median(t))
write_tsv(epdc.medians, "../results/eqtl_dynamic/linear_dQTL/pseudobulk/epdcbin/15bin_medians.tsv")


counts.cmbin <- counts[cm.bins$cell,] %>% 
  aggregate.Matrix(cm.bins$bin.ind, fun="sum") %>% t %>% 
  as_tibble(rownames="gene")
cmbin.drop <- counts.cmbin %>% select(!gene) %>% colSums %>% 
  as_tibble(rownames="sample") %>% filter(value <= 5000) %>% .$sample
counts.cmbin <- counts.cmbin %>% select(-c(cmbin.drop))

counts.epdcbin <- counts[epdc.bins$cell,] %>% 
  aggregate.Matrix(epdc.bins$bin.ind, fun="sum") %>% t %>% 
  as_tibble(rownames="gene")
epdcbin.drop <- counts.epdcbin %>% select(!gene) %>% colSums %>% 
  as_tibble(rownames="sample") %>% filter(value <= 5000) %>% .$sample
counts.epdcbin <- counts.epdcbin %>% select(-c(epdcbin.drop))

# ensure the column names are consistent for the day samples
match_cols <- colnames(bulk.cpm)[colnames(bulk.cpm) %in% colnames(counts.day)]
counts.day <- counts.day %>% relocate(all_of(match_cols))

write_tsv(counts.day, "../data/pseudobulk_counts.day.full.tsv")
write_tsv(counts.type, "../data/pseudobulk_counts.type.full.tsv")
write_tsv(counts.col, "../data/pseudobulk_counts.col.full.tsv")
write_tsv(counts.cmbin, "../data/pseudobulk_counts.cmbin.full.tsv")
write_tsv(counts.epdcbin, "../data/pseudobulk_counts.epdcbin.full.tsv")

# cpm 
cpm.day <- counts.day %>% select(!gene) %>% DGEList %>% calcNormFactors %>% 
  cpm %>% as_tibble %>% mutate(gene=counts.day$gene, .before=1) %>%
  arrange(gene)
cpm.col <- counts.col %>% select(!gene) %>% DGEList %>% calcNormFactors %>% 
  cpm %>% as_tibble %>% mutate(gene=counts.col$gene, .before=1) %>%
  arrange(gene)
cpm.type <- counts.type %>% select(!gene) %>% DGEList %>% calcNormFactors %>% 
  cpm %>% as_tibble %>% mutate(gene=counts.type$gene, .before=1) %>%
  arrange(gene)
cpm.cmbin <- counts.cmbin %>% select(!gene) %>% DGEList %>% calcNormFactors %>% 
  cpm %>% as_tibble %>% mutate(gene=counts.cmbin$gene, .before=1) %>%
  arrange(gene)
cpm.epdcbin <- counts.epdcbin %>% select(!gene) %>% DGEList %>% calcNormFactors %>% 
  cpm %>% as_tibble %>% mutate(gene=counts.epdcbin$gene, .before=1) %>%
  arrange(gene) 

write_tsv(cpm.day, "../data/pseudobulk_cpm.day.full.tsv")
write_tsv(cpm.col, "../data/pseudobulk_cpm.col.full.tsv")
write_tsv(cpm.type, "../data/pseudobulk_cpm.type.full.tsv")
write_tsv(cpm.cmbin, "../data/pseudobulk_cpm.cmbin.full.tsv")
write_tsv(cpm.epdcbin, "../data/pseudobulk_cpm.epdcbin.full.tsv")

# logcpm 
logcpm.day <- counts.day %>% select(!gene) %>% DGEList %>% calcNormFactors %>% 
  cpm(log=T) %>% as_tibble %>% mutate(gene=counts.day$gene) %>% relocate(match_cols) %>%
  arrange(gene)
logcpm.col <- counts.col %>% select(!gene) %>% DGEList %>% calcNormFactors %>% 
  cpm(log=T) %>% as_tibble %>% mutate(gene=counts.col$gene) %>% relocate(gene) %>%
  arrange(gene)
logcpm.type <- counts.type %>% select(!gene) %>% DGEList %>% calcNormFactors %>% 
  cpm(log=T) %>% as_tibble %>% mutate(gene=counts.type$gene) %>% relocate(gene) %>%
  arrange(gene)
logcpm.cmbin <- counts.cmbin %>% select(!gene) %>% DGEList %>% calcNormFactors %>% 
  cpm(log=T) %>% as_tibble %>% mutate(gene=counts.cmbin$gene) %>% relocate(gene) %>%
  arrange(gene)
logcpm.epdcbin <- counts.epdcbin %>% select(!gene) %>% DGEList %>% calcNormFactors %>% 
  cpm(log=T) %>% as_tibble %>% mutate(gene=counts.epdcbin$gene) %>% relocate(gene) %>%
  arrange(gene)


write_tsv(logcpm.day, "../data/pseudobulk_logcpm.day.full.tsv")
write_tsv(logcpm.col, "../data/pseudobulk_logcpm.col.full.tsv")
write_tsv(logcpm.type, "../data/pseudobulk_logcpm.type.full.tsv")
write_tsv(logcpm.cmbin, "../data/pseudobulk_logcpm.cmbin.full.tsv")
write_tsv(logcpm.epdcbin, "../data/pseudobulk_logcpm.epdcbin.full.tsv")

# subset to overlap genes between bulk and pseudobulk
overlap.genes <- intersect(bulk.cpm$gene, cpm.day$gene)
bulk.raw.overlap <- bulk.raw %>% filter(gene %in% overlap.genes)
bulk.cpm.overlap <- bulk.cpm %>% filter(gene %in% overlap.genes)
bulk.rpkm.overlap <- bulk.rpkm %>% filter(gene %in% overlap.genes)
bulk.tpm.overlap <- bulk.tpm %>% filter(gene %in% overlap.genes)
bulk.logtpm.overlap <- bulk.logtpm %>% filter(gene %in% overlap.genes)
counts.day.overlap <- counts.day %>% filter(gene %in% overlap.genes)
counts.col.overlap <- counts.col %>% filter(gene %in% overlap.genes)
counts.type.overlap <- counts.type %>% filter(gene %in% overlap.genes)
cpm.day.overlap <- cpm.day %>% filter(gene %in% overlap.genes)
cpm.col.overlap <- cpm.col %>% filter(gene %in% overlap.genes)
cpm.type.overlap <- cpm.type %>% filter(gene %in% overlap.genes)
logcpm.day.overlap <- logcpm.day %>% filter(gene %in% overlap.genes)
logcpm.col.overlap <- logcpm.col %>% filter(gene %in% overlap.genes)
logcpm.type.overlap <- logcpm.type %>% filter(gene %in% overlap.genes)

write_tsv(bulk.raw.overlap, "../data/bulk_counts.overlap.tsv")
write_tsv(bulk.cpm.overlap, "../data/bulk_cpm.overlap.tsv")
write_tsv(bulk.rpkm.overlap, "../data/bulk_rpkm.overlap.tsv")
write_tsv(bulk.tpm.overlap, "../data/bulk_tpm.overlap.tsv")
write_tsv(bulk.logtpm.overlap, "../data/bulk_logtpm.overlap.tsv")
write_tsv(counts.day.overlap, "../data/pseudobulk_counts.day.overlap.tsv")
write_tsv(counts.col.overlap, "../data/pseudobulk_counts.col.overlap.tsv")
write_tsv(counts.type.overlap, "../data/pseudobulk_counts.type.overlap.tsv")
write_tsv(cpm.day.overlap, "../data/pseudobulk_cpm.day.overlap.tsv")
write_tsv(cpm.col.overlap, "../data/pseudobulk_cpm.col.overlap.tsv")
write_tsv(cpm.type.overlap, "../data/pseudobulk_cpm.type.overlap.tsv")
write_tsv(logcpm.day.overlap, "../data/pseudobulk_logcpm.day.overlap.tsv")
write_tsv(logcpm.col.overlap, "../data/pseudobulk_logcpm.col.overlap.tsv")
write_tsv(logcpm.type.overlap, "../data/pseudobulk_logcpm.type.overlap.tsv")

# run cell line pca on both bulk and pseudobulk, save up to 10 cell line pcs
pseudobulk.clpcs <- cell.line.pca(logcpm.day, npc=10)$cell.line.pcs
pseudobulk.cmbin.clpcs <- cell.line.pca(logcpm.cmbin, npc=10)$cell.line.pcs
pseudobulk.epdcbin.clpcs <- cell.line.pca(logcpm.epdcbin, npc=10)$cell.line.pcs
bulk.clpcs <- cell.line.pca(bulk.logtpm, npc=10)$cell.line.pcs
write_tsv(pseudobulk.clpcs, "../data/dynamic/covariates/pseudobulk.day.full.clpcs.tsv")
write_tsv(pseudobulk.cmbin.clpcs, "../data/dynamic/covariates/pseudobulk.cmbin.full.clpcs.tsv")
write_tsv(pseudobulk.epdcbin.clpcs, "../data/dynamic/covariates/pseudobulk.epdcbin.full.clpcs.tsv")
write_tsv(bulk.clpcs, "../data/dynamic/covariates/bulk.full.clpcs.tsv")

# run regular pca on both bulk and pseudobulk, save up to 50 regular pcs
# this will be used for ieQTL calling, not non-dynamic
pseudobulk.pcs <- regular.pca(logcpm.col, 50)$u
pseudobulk.cmbin.pcs <- regular.pca(logcpm.cmbin, 50)$u
pseudobulk.epdcbin.pcs <- regular.pca(logcpm.epdcbin, 50)$u
pseudobulk.day.pcs <- regular.pca(logcpm.day, 50)$u
bulk.pcs <- regular.pca(bulk.logtpm, 50)$u
write_tsv(pseudobulk.pcs, "../data/dynamic/covariates/pseudobulk.col.full.pcs.tsv")
write_tsv(pseudobulk.cmbin.pcs, "../data/dynamic/covariates/pseudobulk.cmbin.full.pcs.tsv")
write_tsv(pseudobulk.epdcbin.pcs, "../data/dynamic/covariates/pseudobulk.epdcbin.full.pcs.tsv")
write_tsv(pseudobulk.day.pcs, "../data/dynamic/covariates/pseudobulk.day.full.pcs.tsv")
write_tsv(bulk.pcs, "../data/dynamic/covariates/bulk.full.pcs.tsv")

# subset to each day, and each cell type, and save up to 5 pcs
# this will be used for non-dynamic eQTL calling
source("../code/cell_line_pca.R")
logcpm.day <- logcpm.day %>%
  column_to_rownames("gene") %>% 
  t %>% as_tibble(rownames="sample") %>%
  mutate(day=str_extract(sample, "[^_]+$"))
for (d in unique(logcpm.day$day)) {
  expr <- logcpm.day %>% 
    filter(day==d) %>%
    mutate(sample=str_sub(sample, 1, 5)) %>%
    select(!day) %>%
    write_tsv(paste0("../data/static/pseudobulk/day/day", d, "/expression.tsv"))
  cvrt <- expr %>%
    column_to_rownames("sample") %>%
    t %>% as_tibble(rownames="gene") %>%
    regular.pca %>%
    .$u %>%
    write_tsv(paste0("../data/static/pseudobulk/day/day", d, "/covariates.tsv"))
}

logcpm.type <- logcpm.type %>% 
  column_to_rownames("gene") %>% 
  t %>% as_tibble(rownames="sample") %>%
  mutate(type=str_extract(sample, "[^_]+$"))
for (t in unique(logcpm.type$type)) {
  expr <- logcpm.type %>% 
    filter(type==t) %>%
    mutate(sample=str_sub(sample, 1, 5)) %>%
    select(!type) %>%
    write_tsv(paste0("../data/static/pseudobulk/type/", t, "/expression.tsv"))
  cvrt <- expr %>%
    column_to_rownames("sample") %>%
    t %>% as_tibble(rownames="gene") %>%
    regular.pca %>%
    .$u %>%
    write_tsv(paste0("../data/static/pseudobulk/type/", t, "/covariates.tsv"))
}

logcpm.cmbin <- logcpm.cmbin %>% 
  column_to_rownames("gene") %>% 
  t %>% as_tibble(rownames="sample") %>%
  mutate(bin=str_extract(sample, "[^_]+$"))
for (b in unique(logcpm.cmbin$bin)) {
  expr <- logcpm.cmbin %>% 
    filter(bin==b) %>%
    mutate(sample=str_sub(sample, 1, 5)) %>%
    select(!bin) %>%
    write_tsv(paste0("../data/static/pseudobulk/cmbin/bin", b, "/expression.tsv"))
  cvrt <- expr %>%
    column_to_rownames("sample") %>%
    t %>% as_tibble(rownames="gene") %>%
    regular.pca %>%
    .$u %>%
    write_tsv(paste0("../data/static/pseudobulk/cmbin/bin", b, "/covariates.tsv"))
}

bulk.logtpm <- bulk.logtpm %>%
  column_to_rownames("gene") %>% 
  t %>% as_tibble(rownames="sample") %>%
  mutate(day=str_extract(sample, "[^_]+$"))
for (d in unique(bulk.logtpm$day)) {
  expr <- bulk.logtpm %>% 
    filter(day==d) %>%
    mutate(sample=str_sub(sample, 1, 5)) %>%
    select(!day) %>%
    write_tsv(paste0("../data/static/bulk/day/day", d, "/expression.tsv"))
  cvrt <- expr %>%
    column_to_rownames("sample") %>%
    t %>% as_tibble(rownames="gene") %>%
    regular.pca %>%
    .$u %>%
    write_tsv(paste0("../data/static/bulk/day/day", d, "/covariates.tsv"))
}

