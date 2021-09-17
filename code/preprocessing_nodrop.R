library(assertthat)
library(tidyverse)
library(edgeR)
library(Seurat)
library(Matrix)
library(Matrix.utils)
library(corrplot)
set.seed(2021)
source("/project2/gilad/jpopp/sc-dynamic-eqtl/code/cell_line_pca.R")

# this file is for examining the impact of dropping samples on dynamic eQTL calling

# HELPER FUNCTIONS
filter_by_depth <- function(c, libsize.cutoff=1e5) {
  c <- select(c, !gene) 
  drop_samples <- DGEList(counts=c) %>%
    .$sample %>% as_tibble(rownames="sample") %>% 
    filter(lib.size<libsize.cutoff) %>% select(sample)
  drop_samples
}

counts2cpm <- function(c, lognorm=F) {
  g <- c$gene
  c <- select(c, !gene)
  
  cpm <- DGEList(counts=c) %>% calcNormFactors(method="TMMwsp") %>% cpm(log=lognorm) %>%
    as_tibble %>% mutate(gene=g, .before=1) %>% arrange(gene)
  cpm
}

ensg2hgnc <- function(cpm, trans=translator) {
  hgnc <- cpm %>% inner_join(trans, by="gene") %>% select(!gene) %>% 
    rename(gene=gene_name) %>% relocate(gene) %>% arrange(gene)
  hgnc
}

# subset to cardiomyocytes (CM)
sc_cm <- readRDS("../data/seurat.cm.rds")
counts_cm <- sc_cm@assays[["SCT"]]@counts %>% t

# bin CM by pseudotime
nbins <- c(16)
for (nb in nbins) {
  cm.bins <- as_tibble(sc_cm$pseudotime, rownames="cell") %>%
    rename(t=value) %>% arrange(t) %>%
    rowid_to_column("time_order") %>% 
    mutate(bin=floor(nb*time_order/nrow(.)), .keep="unused") %>%
    mutate(bin=sapply(bin, function(x){min(x, nb-1)}))
  cm.binind <- as_tibble(sc_cm$individual, rownames="cell") %>%
    rename(ind=value) %>%
    left_join(cm.bins, by="cell") %>%
    mutate(binind=paste(ind, bin, sep="_")) 
  cm.ncells <- cm.binind %>%
    group_by(binind) %>%
    count %>%
    write_tsv(paste0("../results/eqtl_dynamic/linear_dQTL/pseudobulk-cm-nodrop/bin", nb, "/bin_ncells.tsv"))
  cm.medians <- cm.binind %>%
    group_by(binind) %>%
    summarize(t=median(t)) %>%
    write_tsv(paste0("../results/eqtl_dynamic/linear_dQTL/pseudobulk-cm-nodrop/bin", nb, "/bin_medians.tsv"))
  cm.counts.bin <- counts_cm %>% aggregate.Matrix(cm.binind$binind, fun="sum") %>% t %>%
    as_tibble(rownames="gene") %>% arrange(gene) %>% write_tsv(paste0("../data/pseudobulk-cm-nodrop/bin", nb, "/counts.tsv"))
  cm.libsize <- cm.counts.bin %>% 
    select(-c(gene)) %>%
    DGEList %>%
    calcNormFactors(method="TMMwsp") %>%
    .$samples %>%
    as_tibble(rownames="sample") %>%
    select(!group) %>%
    write_tsv(paste0("../results/eqtl_dynamic/linear_dQTL/pseudobulk-cm-nodrop/bin", nb, "/bin_libsize.tsv"))
  cm.cpm.bin <- cm.counts.bin %>% counts2cpm %>% write_tsv(paste0("../data/pseudobulk-cm-nodrop/bin", nb, "/cpm.tsv"))
  cm.logcpm.bin <- cm.counts.bin %>% counts2cpm(lognorm=T) %>% write_tsv(paste0("../data/pseudobulk-cm-nodrop/bin", nb, "/logcpm.tsv"))
  
  pseudobulk.cm.bin.clpcs <- cell.line.pca(cm.logcpm.bin, npc=10)$cell.line.pcs %>% write_tsv(paste0("../data/pseudobulk-cm-nodrop/bin", nb, "/clpcs.tsv"))
  pseudobulk.cm.bin.pcs <- regular.pca(cm.logcpm.bin, 50)$u %>% write_tsv(paste0("../data/pseudobulk-cm-nodrop/bin", nb, "/pcs.tsv"))
}

# subset to cardiomyocytes (CF)
sc_cf <- readRDS("../data/seurat.cf.rds")
counts_cf <- sc_cf@assays[["SCT"]]@counts %>% t

# bin CF by pseudotime
nbins <- c(16)
for (nb in nbins) {
  cf.bins <- as_tibble(sc_cf$pseudotime, rownames="cell") %>%
    rename(t=value) %>% arrange(t) %>%
    rowid_to_column("time_order") %>% 
    mutate(bin=floor(nb*time_order/nrow(.)), .keep="unused") %>%
    mutate(bin=sapply(bin, function(x){min(x, nb-1)}))
  cf.binind <- as_tibble(sc_cf$individual, rownames="cell") %>%
    rename(ind=value) %>%
    left_join(cf.bins, by="cell") %>%
    mutate(binind=paste(ind, bin, sep="_")) 
  cf.ncells <- cf.binind %>%
    group_by(binind) %>%
    count %>%
    write_tsv(paste0("../results/eqtl_dynamic/linear_dQTL/pseudobulk-cf-nodrop/bin", nb, "/bin_ncells.tsv"))
  cf.medians <- cf.binind %>%
    group_by(binind) %>%
    summarize(t=median(t)) %>%
    write_tsv(paste0("../results/eqtl_dynamic/linear_dQTL/pseudobulk-cf-nodrop/bin", nb, "/bin_medians.tsv"))
  cf.counts.bin <- counts_cf %>% aggregate.Matrix(cf.binind$binind, fun="sum") %>% t %>%
    as_tibble(rownames="gene") %>% arrange(gene) %>% write_tsv(paste0("../data/pseudobulk-cf-nodrop/bin", nb, "/counts.tsv"))
  cf.libsize <- cf.counts.bin %>% 
    select(-c(gene)) %>%
    DGEList %>%
    calcNormFactors(method="TMMwsp") %>%
    .$samples %>%
    as_tibble(rownames="sample") %>%
    select(!group) %>%
    write_tsv(paste0("../results/eqtl_dynamic/linear_dQTL/pseudobulk-cf-nodrop/bin", nb, "/bin_libsize.tsv"))
  cf.cpm.bin <- cf.counts.bin %>% counts2cpm %>% write_tsv(paste0("../data/pseudobulk-cf-nodrop/bin", nb, "/cpm.tsv"))
  cf.logcpm.bin <- cf.counts.bin %>% counts2cpm(lognorm=T) %>% write_tsv(paste0("../data/pseudobulk-cf-nodrop/bin", nb, "/logcpm.tsv"))
  
  pseudobulk.cf.bin.clpcs <- cell.line.pca(cf.logcpm.bin, npc=10)$cell.line.pcs %>% write_tsv(paste0("../data/pseudobulk-cf-nodrop/bin", nb, "/clpcs.tsv"))
  pseudobulk.cf.bin.pcs <- regular.pca(cf.logcpm.bin, 50)$u %>% write_tsv(paste0("../data/pseudobulk-cf-nodrop/bin", nb, "/pcs.tsv"))
}