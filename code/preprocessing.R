library(assertthat)
library(tidyverse)
library(edgeR)
library(Seurat)
library(Matrix)
library(Matrix.utils)
library(corrplot)
set.seed(2021)
source("/project2/gilad/jpopp/sc-dynamic-eqtl/code/cell_line_pca.R")

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

# BULK 
# load the bulk data, convert ENSG to HGNC, convert to RPKM (to compare to Ben's) and TPM
translator <- read_tsv("/project2/gilad/jpopp/cm_eqtl/data/gencode.v26.annotation.gene.txt") %>%
  mutate(gene=str_extract(gene_id, "[^.]+")) %>%
  filter(!duplicated(gene)) %>% select(gene, gene_name)
bulk <- read_delim("/project2/gilad/bstrober/ipsc_differentiation_19_lines/preprocess/processed_total_expression/raw_counts.txt",
                   delim=" ") %>% rename(gene=`Gene_id`)
gene_len <- read_tsv("/project2/gilad/bstrober/ipsc_differentiation_19_lines/preprocess/processed_total_expression/gene_lengths.txt") %>%
  rename(gene=GeneID)
counts <- bulk %>% inner_join(gene_len, by="gene") %>% ensg2hgnc %>% select(!Length)
lengths <- bulk %>% inner_join(gene_len, by="gene") %>% ensg2hgnc %>% select(Length)
rpk <- as_tibble(select(counts, !gene)/lengths$Length) %>% mutate(gene=counts$gene)

drop.samples <- counts %>% filter_by_depth %>% write_tsv("../data/bulk/day/dropped_samples.tsv")
bulk.counts <- counts %>% select(setdiff(colnames(.), drop.samples$sample)) %>% write_tsv("../data/bulk/day/counts.tsv")
bulk.rpkm <- DGEList(counts=select(counts, !gene), genes=lengths) %>% calcNormFactors(method="TMMwsp") %>% rpkm %>%
  as_tibble %>% mutate(gene=counts$gene, .before=1) %>% write_tsv("../data/bulk/day/rpkm.tsv")
bulk.cpm <- counts2cpm(counts) %>% write_tsv("../data/bulk/day/cpm.tsv")
bulk.tpm <- counts2cpm(rpk) %>% write_tsv("../data/bulk/day/tpm.tsv")
bulk.logtpm <- counts2cpm(rpk, lognorm=T) %>% write_tsv("../data/bulk/day/logtpm.tsv")
bulk.clpcs <- cell.line.pca(bulk.logtpm, npc=10)$cell.line.pcs %>% write_tsv("../data/bulk/day/clpcs.tsv")
bulk.pcs <- regular.pca(bulk.logtpm, 50)$u %>% write_tsv("../data/bulk/day/pcs.tsv")

# subset to the 7 days that are re-collected in the single cell experiment
pb.days <- c(0, 1, 3, 5, 7, 11, 15)
counts7 <- counts %>% select(c(gene, ends_with(paste0("_", pb.days)))) 
rpk7 <- as_tibble(select(counts7, !gene)/lengths$Length) %>% mutate(gene=counts7$gene)

drop.samples <- counts7 %>% filter_by_depth %>% write_tsv("../data/bulk7/day/dropped_samples.tsv")
bulk7.counts <- counts7 %>% select(setdiff(colnames(.), drop.samples$sample)) %>% write_tsv("../data/bulk7/day/counts.tsv")
bulk7.rpkm <- DGEList(counts=select(counts7, !gene), genes=lengths) %>% calcNormFactors %>% rpkm %>%
  as_tibble %>% mutate(gene=counts7$gene, .before=1) %>% write_tsv("../data/bulk7/day/rpkm.tsv")
bulk7.cpm <- counts2cpm(counts7) %>% write_tsv("../data/bulk7/day/cpm.tsv")
bulk7.tpm <- counts2cpm(rpk7) %>% write_tsv("../data/bulk7/day/tpm.tsv")
bulk7.logtpm <- counts2cpm(rpk7, lognorm=T) %>% write_tsv("../data/bulk7/day/logtpm.tsv")
bulk7.clpcs <- cell.line.pca(bulk7.logtpm, npc=10)$cell.line.pcs %>% write_tsv("../data/bulk7/day/clpcs.tsv")
bulk7.pcs <- regular.pca(bulk7.logtpm, 50)$u %>% write_tsv("../data/bulk7/day/pcs.tsv")

# subset each of these to 
# PSEUDOBULK
sc <- readRDS("../data/seurat.annotated.rds")
counts <- sc@assays[["SCT"]]@counts %>% t

# aggregation
# bin all cells by type
type.ind <- paste(sc$individual, sc$type, sep="_")
counts.type <- counts %>% aggregate.Matrix(type.ind, fun="sum") %>% t %>% 
  as_tibble(rownames="gene") %>% arrange(gene)
drop.samples <- counts.type %>% filter_by_depth %>% write_tsv("../data/pseudobulk/type/dropped_samples.tsv")
counts.type <- counts.type %>% select(setdiff(colnames(.), drop.samples$sample)) %>% write_tsv("../data/pseudobulk/type/counts.tsv")
cpm.type <- counts.type %>% counts2cpm %>% write_tsv("../data/pseudobulk/type/cpm.tsv")
logcpm.type <- counts.type %>% counts2cpm(lognorm=T) %>% write_tsv("../data/pseudobulk/type/logcpm.tsv")

# bin all cells by collection day
day.ind <- paste(sc$individual, sc$diffday, sep="_") %>% str_replace("day", "")
counts.day <- counts %>% aggregate.Matrix(day.ind, fun="sum") %>% t %>%
  as_tibble(rownames="gene") %>% arrange(gene)
drop.samples <- counts.day %>% filter_by_depth %>% write_tsv("../data/pseudobulk/day/dropped_samples.tsv")
counts.day <- counts.day %>% select(setdiff(colnames(.), drop.samples$sample)) %>% write_tsv("../data/pseudobulk/day/counts.tsv")
cpm.day <- counts.day %>% counts2cpm %>% write_tsv("../data/pseudobulk/day/cpm.tsv")
logcpm.day <- counts.day %>% counts2cpm(lognorm=T) %>% write_tsv("../data/pseudobulk/day/logcpm.tsv")
pseudobulk.day.clpcs <- cell.line.pca(logcpm.day, npc=10)$cell.line.pcs %>% write_tsv("../data/pseudobulk/day/clpcs.tsv")
pseudobulk.day.pcs <- regular.pca(logcpm.day, 50)$u %>% write_tsv("../data/pseudobulk/day/pcs.tsv")

# bin all cells by collection
col.ind <- paste(sc$individual, sc$orig.ident, sep="_")
counts.col <- counts %>% aggregate.Matrix(col.ind, fun="sum") %>% t %>% 
  as_tibble(rownames="gene") %>% arrange(gene)
drop.samples <- counts.col %>% filter_by_depth %>% write_tsv("../data/pseudobulk/col/dropped_samples.tsv")
counts.col <- counts.col %>% select(setdiff(colnames(.), drop.samples$sample)) %>% write_tsv("../data/pseudobulk/col/counts.tsv")
cpm.col <- counts.col %>% counts2cpm %>% write_tsv("../data/pseudobulk/col/cpm.tsv")
logcpm.col <- counts.col %>% counts2cpm(lognorm=T) %>% write_tsv("../data/pseudobulk/col/logcpm.tsv")
pseudobulk.col.pcs <- regular.pca(logcpm.col, 50)$u %>% write_tsv("../data/pseudobulk/col/pcs.tsv")

# subset to cardiomyocytes (CM)
sc_cm <- readRDS("../data/seurat.cm.rds")
counts_cm <- sc_cm@assays[["SCT"]]@counts %>% t

# bin CM by collection day
cm.day.ind <- paste(sc_cm$individual, sc_cm$diffday, sep="_") %>% str_replace("day", "")
cm.counts.day <- counts_cm %>% aggregate.Matrix(cm.day.ind, fun="sum") %>% t %>%
  as_tibble(rownames="gene") %>% arrange(gene)
drop.samples <- cm.counts.day %>% filter_by_depth %>% write_tsv("../data/pseudobulk-cm/day/dropped_samples.tsv")
cm.counts.day <- cm.counts.day %>% select(setdiff(colnames(.), drop.samples$sample)) %>% write_tsv("../data/pseudobulk-cm/day/counts.tsv")
cm.cpm.day <- cm.counts.day %>% counts2cpm %>% write_tsv("../data/pseudobulk-cm/day/cpm.tsv")
cm.logcpm.day <- cm.counts.day %>% counts2cpm(lognorm=T) %>% write_tsv("../data/pseudobulk-cm/day/logcpm.tsv")
pseudobulk.cm.day.clpcs <- cell.line.pca(cm.logcpm.day, npc=10)$cell.line.pcs %>% write_tsv("../data/pseudobulk-cm/day/clpcs.tsv")
pseudobulk.cm.day.pcs <- regular.pca(cm.logcpm.day, 50)$u %>% write_tsv("../data/pseudobulk-cm/day/pcs.tsv")

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
    write_tsv(paste0("../results/eqtl_dynamic/linear_dQTL/pseudobulk-cm/bin", nb, "/bin_ncells.tsv"))
  cm.medians <- cm.binind %>%
    group_by(binind) %>%
    summarize(t=median(t)) %>%
    write_tsv(paste0("../results/eqtl_dynamic/linear_dQTL/pseudobulk-cm/bin", nb, "/bin_medians.tsv"))
  cm.counts.bin <- counts_cm %>% aggregate.Matrix(cm.binind$binind, fun="sum") %>% t %>%
    as_tibble(rownames="gene") %>% arrange(gene)
  drop.samples <- cm.counts.bin %>% filter_by_depth %>% write_tsv(paste0("../data/pseudobulk-cm/bin", nb, "/dropped_samples.tsv"))
  cm.counts.bin <- cm.counts.bin %>% select(setdiff(colnames(.), drop.samples$sample)) %>% write_tsv(paste0("../data/pseudobulk-cm/bin", nb, "/counts.tsv"))
  cm.libsize <- cm.counts.bin %>% 
    select(setdiff(colnames(.), drop.samples$sample)) %>%
    select(-c(gene)) %>%
    DGEList %>%
    calcNormFactors(method="TMMwsp") %>%
    .$samples %>%
    as_tibble(rownames="sample") %>%
    select(!group) %>%
    write_tsv(paste0("../results/eqtl_dynamic/linear_dQTL/pseudobulk-cm/bin", nb, "/bin_libsize.tsv"))
  cm.cpm.bin <- cm.counts.bin %>% counts2cpm %>% write_tsv(paste0("../data/pseudobulk-cm/bin", nb, "/cpm.tsv"))
  cm.logcpm.bin <- cm.counts.bin %>% counts2cpm(lognorm=T) %>% write_tsv(paste0("../data/pseudobulk-cm/bin", nb, "/logcpm.tsv"))
  
  pseudobulk.cm.bin.clpcs <- cell.line.pca(cm.logcpm.bin, npc=10)$cell.line.pcs %>% write_tsv(paste0("../data/pseudobulk-cm/bin", nb, "/clpcs.tsv"))
  pseudobulk.cm.bin.pcs <- regular.pca(cm.logcpm.bin, 50)$u %>% write_tsv(paste0("../data/pseudobulk-cm/bin", nb, "/pcs.tsv"))
}

# subset to cardiomyocytes (CF)
sc_cf <- readRDS("../data/seurat.cf.rds")
counts_cf <- sc_cf@assays[["SCT"]]@counts %>% t

# bin CF by collection day
cf.day.ind <- paste(sc_cf$individual, sc_cf$diffday, sep="_") %>% str_replace("day", "")
cf.counts.day <- counts_cf %>% aggregate.Matrix(cf.day.ind, fun="sum") %>% t %>%
  as_tibble(rownames="gene") %>% arrange(gene) 
drop.samples <- cf.counts.day %>% filter_by_depth %>% write_tsv("../data/pseudobulk-cf/day/dropped_samples.tsv")
cf.counts.day <- cf.counts.day %>% select(setdiff(colnames(.), drop.samples$sample)) %>% write_tsv("../data/pseudobulk-cf/day/counts.tsv")
cf.cpm.day <- cf.counts.day %>% counts2cpm %>% write_tsv("../data/pseudobulk-cf/day/cpm.tsv")
cf.logcpm.day <- cf.counts.day %>% counts2cpm(lognorm=T) %>% write_tsv("../data/pseudobulk-cf/day/logcpm.tsv")
pseudobulk.cf.day.clpcs <- cell.line.pca(cf.logcpm.day, npc=10)$cell.line.pcs %>% write_tsv("../data/pseudobulk-cf/day/clpcs.tsv")
pseudobulk.cf.day.pcs <- regular.pca(cf.logcpm.day, 50)$u %>% write_tsv("../data/pseudobulk-cf/day/pcs.tsv")

# bin CF by pseudotime
nbins <- c(7, 16, 25, 30)
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
    write_tsv(paste0("../results/eqtl_dynamic/linear_dQTL/pseudobulk-cf/bin", nb, "/bin_ncells.tsv"))
  cf.medians <- cf.binind %>%
    group_by(binind) %>%
    summarize(t=median(t)) %>%
    write_tsv(paste0("../results/eqtl_dynamic/linear_dQTL/pseudobulk-cf/bin", nb, "/bin_medians.tsv"))
  cf.counts.bin <- counts_cf %>% aggregate.Matrix(cf.binind$binind, fun="sum") %>% t %>%
    as_tibble(rownames="gene") %>% arrange(gene) 
  drop.samples <- cf.counts.bin %>% filter_by_depth %>% write_tsv(paste0("../data/pseudobulk-cf/bin", nb, "/dropped_samples.tsv"))
  cf.counts.bin <- cf.counts.bin %>% select(setdiff(colnames(.), drop.samples$sample)) %>% write_tsv(paste0("../data/pseudobulk-cf/bin", nb, "/counts.tsv"))
  cf.libsize <- cf.counts.bin %>% 
    select(setdiff(colnames(.), drop.samples$sample)) %>%
    select(-c(gene)) %>%
    DGEList %>%
    calcNormFactors(method="TMMwsp") %>%
    .$samples %>%
    as_tibble(rownames="sample") %>%
    select(!group) %>%
    write_tsv(paste0("../results/eqtl_dynamic/linear_dQTL/pseudobulk-cf/bin", nb, "/bin_libsize.tsv"))
  cf.cpm.bin <- cf.counts.bin %>% counts2cpm %>% write_tsv(paste0("../data/pseudobulk-cf/bin", nb, "/cpm.tsv"))
  cf.logcpm.bin <- cf.counts.bin %>% counts2cpm(lognorm=T) %>% write_tsv(paste0("../data/pseudobulk-cf/bin", nb, "/logcpm.tsv"))
  
  pseudobulk.cf.bin.clpcs <- cell.line.pca(cf.logcpm.bin, npc=10)$cell.line.pcs %>% write_tsv(paste0("../data/pseudobulk-cf/bin", nb, "/clpcs.tsv"))
  pseudobulk.cf.bin.pcs <- regular.pca(cf.logcpm.bin, 50)$u %>% write_tsv(paste0("../data/pseudobulk-cf/bin", nb, "/pcs.tsv"))
}
# 
# subset to overlap genes between bulk and pseudobulk
# overlap.genes <- intersect(bulk.cpm$gene, cpm.day$gene)
# bulk.raw.overlap <- bulk.raw %>% filter(gene %in% overlap.genes)
# bulk.cpm.overlap <- bulk.cpm %>% filter(gene %in% overlap.genes)
# bulk.rpkm.overlap <- bulk.rpkm %>% filter(gene %in% overlap.genes)
# bulk.tpm.overlap <- bulk.tpm %>% filter(gene %in% overlap.genes)
# bulk.logtpm.overlap <- bulk.logtpm %>% filter(gene %in% overlap.genes)
# counts.day.overlap <- counts.day %>% filter(gene %in% overlap.genes)
# counts.col.overlap <- counts.col %>% filter(gene %in% overlap.genes)
# counts.type.overlap <- counts.type %>% filter(gene %in% overlap.genes)
# cpm.day.overlap <- cpm.day %>% filter(gene %in% overlap.genes)
# cpm.col.overlap <- cpm.col %>% filter(gene %in% overlap.genes)
# cpm.type.overlap <- cpm.type %>% filter(gene %in% overlap.genes)
# logcpm.day.overlap <- logcpm.day %>% filter(gene %in% overlap.genes)
# logcpm.col.overlap <- logcpm.col %>% filter(gene %in% overlap.genes)
# logcpm.type.overlap <- logcpm.type %>% filter(gene %in% overlap.genes)
# 
# write_tsv(bulk.raw.overlap, "../data/bulk_counts.overlap.tsv")
# write_tsv(bulk.cpm.overlap, "../data/bulk_cpm.overlap.tsv")
# write_tsv(bulk.rpkm.overlap, "../data/bulk_rpkm.overlap.tsv")
# write_tsv(bulk.tpm.overlap, "../data/bulk_tpm.overlap.tsv")
# write_tsv(bulk.logtpm.overlap, "../data/bulk_logtpm.overlap.tsv")
# write_tsv(counts.day.overlap, "../data/pseudobulk_counts.day.overlap.tsv")
# write_tsv(counts.col.overlap, "../data/pseudobulk_counts.col.overlap.tsv")
# write_tsv(counts.type.overlap, "../data/pseudobulk_counts.type.overlap.tsv")
# write_tsv(cpm.day.overlap, "../data/pseudobulk_cpm.day.overlap.tsv")
# write_tsv(cpm.col.overlap, "../data/pseudobulk_cpm.col.overlap.tsv")
# write_tsv(cpm.type.overlap, "../data/pseudobulk_cpm.type.overlap.tsv")
# write_tsv(logcpm.day.overlap, "../data/pseudobulk_logcpm.day.overlap.tsv")
# write_tsv(logcpm.col.overlap, "../data/pseudobulk_logcpm.col.overlap.tsv")
# write_tsv(logcpm.type.overlap, "../data/pseudobulk_logcpm.type.overlap.tsv")

# run cell line pca on both bulk and pseudobulk, save up to 10 cell line pcs


# run regular pca on both bulk and pseudobulk, save up to 50 regular pcs
# this will be used for ieQTL calling, not non-dynamic


# subset to each day, and each cell type, and save up to 5 pcs
# this will be used for non-dynamic eQTL calling
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

bulk7.logtpm <- bulk7.logtpm %>%
  column_to_rownames("gene") %>%
  t %>% as_tibble(rownames="sample") %>%
  mutate(day=str_extract(sample, "[^_]+$"))
for (d in unique(bulk7.logtpm$day)) {
  expr <- bulk7.logtpm %>%
    filter(day==d) %>%
    mutate(sample=str_sub(sample, 1, 5)) %>%
    select(!day) %>%
    write_tsv(paste0("../data/static/bulk7/day/day", d, "/expression.tsv"))
  cvrt <- expr %>%
    column_to_rownames("sample") %>%
    t %>% as_tibble(rownames="gene") %>%
    regular.pca %>%
    .$u %>%
    write_tsv(paste0("../data/static/bulk7/day/day", d, "/covariates.tsv"))
}

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

cm.logcpm.bin <- cm.logcpm.bin %>%
  column_to_rownames("gene") %>%
  t %>% as_tibble(rownames="sample") %>%
  mutate(bin=str_extract(sample, "[^_]+$"))
for (b in unique(cm.logcpm.bin$bin)) {
  expr <- cm.logcpm.bin %>%
    filter(bin==b) %>%
    mutate(sample=str_sub(sample, 1, 5)) %>%
    select(!bin) %>%
    write_tsv(paste0("../data/static/pseudobulk-cm/bin15/bin", b, "/expression.tsv"))
  cvrt <- expr %>%
    column_to_rownames("sample") %>%
    t %>% as_tibble(rownames="gene") %>%
    regular.pca %>%
    .$u %>%
    write_tsv(paste0("../data/static/pseudobulk-cm/bin15/bin", b, "/covariates.tsv"))
}

cf.logcpm.bin <- cf.logcpm.bin %>%
  column_to_rownames("gene") %>%
  t %>% as_tibble(rownames="sample") %>%
  mutate(bin=str_extract(sample, "[^_]+$"))
for (b in unique(cf.logcpm.bin$bin)) {
  expr <- cf.logcpm.bin %>%
    filter(bin==b) %>%
    mutate(sample=str_sub(sample, 1, 5)) %>%
    select(!bin) %>%
    write_tsv(paste0("../data/static/pseudobulk-cf/bin15/bin", b, "/expression.tsv"))
  cvrt <- expr %>%
    column_to_rownames("sample") %>%
    t %>% as_tibble(rownames="gene") %>%
    regular.pca %>%
    .$u %>%
    write_tsv(paste0("../data/static/pseudobulk-cf/bin15/bin", b, "/covariates.tsv"))
}

