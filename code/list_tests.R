library(tidyverse)

gene_locs <- read_tsv("../data/gene_locs.filtered.tsv")
snp_locs <- read_tsv("../data/snp_locs.filtered.tsv")

cis.dist <- 50000
get.snps <- function(i) {
  tss <- gene_locs[i, "start_pos"][[1]]
  g.snps <- filter(snp_locs, chr==gene_locs[i, "chr"][[1]] & pos>=(tss-cis.dist) & pos<=(tss+cis.dist))
  d2tss <- abs(g.snps$pos-tss)
  tibble(gene=gene_locs[i, "gene_name"][[1]], snp=g.snps$snp, dist2tss=d2tss)
}

all.tests <- bind_rows(lapply(1:nrow(gene_locs), get.snps))
write_tsv(all.tests, "../data/gv_pairs.filtered.50k.tsv")
