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

bulk.days <- as.character(seq(0, 15))
pseudobulk.days <- c(0, 1, 3, 5, 7, 11, 15)
pseudobulk.types <- c("IPSC", "MES", "CMES", "PROG", "CM", "CF", "UNK")

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

for (d in bulk.days) {
  bulk.counts <- counts %>% select(c(gene, ends_with(paste0("_", d)))) %>% rename_with(function(x){if_else(x=="gene", true=x, false=str_sub(x, 1, 5))}) %>% write_tsv(paste0("../data/static/bulk/day/day", d, "/counts.tsv"))
  bulk.tpm <- counts2cpm(select(rpk, c(gene, ends_with(paste0("_", d))))) %>% rename_with(function(x){if_else(x=="gene", true=x, false=str_sub(x, 1, 5))}) %>% write_tsv(paste0("../data/static/bulk/day/day", d, "/tpm.tsv"))
  bulk.logtpm <- counts2cpm(select(rpk, c(gene, ends_with(paste0("_", d)))), lognorm=T) %>% rename_with(function(x){if_else(x=="gene", true=x, false=str_sub(x, 1, 5))}) %>% write_tsv(paste0("../data/static/bulk/day/day", d, "/logtpm.tsv"))
  exp <- bulk.logtpm %>% column_to_rownames("gene") %>% t %>% as_tibble(rownames="sample") %>% write_tsv(paste0("../data/static/bulk/day/day", d, "/expression.tsv"))
  bulk.pcs <- regular.pca(bulk.logtpm, 5)$u %>% write_tsv(paste0("../data/static/bulk/day/day", d, "/pcs.tsv"))
}

# PSEUDOBULK
sc <- readRDS("../data/seurat.clustered.rds")
counts <- sc@assays[["SCT"]]@counts %>% t

# aggregation by type
type.ind <- paste(sc$individual, sc$type, sep="_")
counts.type <- counts %>% aggregate.Matrix(type.ind, fun="sum") %>% t %>% 
  as_tibble(rownames="gene") %>% arrange(gene)
for (d in pseudobulk.types) {
  type.counts <- counts.type %>% select(c(gene, ends_with(paste0("_", d)))) %>% rename_with(function(x){if_else(x=="gene", true=x, false=str_sub(x, 1, 5))}) %>% write_tsv(paste0("../data/static/pseudobulk/type/", d, "/counts.tsv"))
  type.cpm <- counts.type %>% select(c(gene, ends_with(paste0("_", d)))) %>% rename_with(function(x){if_else(x=="gene", true=x, false=str_sub(x, 1, 5))}) %>% counts2cpm %>% write_tsv(paste0("../data/static/pseudobulk/type/", d, "/cpm.tsv"))
  type.logcpm <- counts.type %>% select(c(gene, ends_with(paste0("_", d)))) %>% rename_with(function(x){if_else(x=="gene", true=x, false=str_sub(x, 1, 5))}) %>% counts2cpm(lognorm=T) %>% write_tsv(paste0("../data/static/pseudobulk/type/", d, "/logcpm.tsv"))
  exp <- type.logcpm %>% column_to_rownames("gene") %>% t %>% as_tibble(rownames="sample") %>% write_tsv(paste0("../data/static/pseudobulk/type/", d, "/expression.tsv"))
  type.pcs <- regular.pca(type.logcpm, 5)$u %>% write_tsv(paste0("../data/static/pseudobulk/type/", d, "/pcs.tsv"))
}

# aggregation by day
day.ind <- paste(sc$individual, sc$diffday, sep="_") %>% str_replace("day", "")
counts.day <- counts %>% aggregate.Matrix(day.ind, fun="sum") %>% t %>%
  as_tibble(rownames="gene") %>% arrange(gene)
for (d in pseudobulk.days) {
  day.counts <- counts.day %>% select(c(gene, ends_with(paste0("_", d)))) %>% rename_with(function(x){if_else(x=="gene", true=x, false=str_sub(x, 1, 5))}) %>% write_tsv(paste0("../data/static/pseudobulk/day/day", d, "/counts.tsv"))
  day.cpm <- counts.day %>% select(c(gene, ends_with(paste0("_", d)))) %>% rename_with(function(x){if_else(x=="gene", true=x, false=str_sub(x, 1, 5))}) %>% counts2cpm %>% write_tsv(paste0("../data/static/pseudobulk/day/day", d, "/cpm.tsv"))
  day.logcpm <- counts.day %>% select(c(gene, ends_with(paste0("_", d)))) %>% rename_with(function(x){if_else(x=="gene", true=x, false=str_sub(x, 1, 5))}) %>% counts2cpm(lognorm=T) %>% write_tsv(paste0("../data/static/pseudobulk/day/day", d, "/logcpm.tsv"))
  exp <- day.logcpm %>% column_to_rownames("gene") %>% t %>% as_tibble(rownames="sample") %>% write_tsv(paste0("../data/static/pseudobulk/day/day", d, "/expression.tsv"))
  day.pcs <- regular.pca(day.logcpm, 5)$u %>% write_tsv(paste0("../data/static/pseudobulk/day/day", d, "/pcs.tsv"))
}

