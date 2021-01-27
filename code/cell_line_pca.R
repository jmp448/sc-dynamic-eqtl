library(tidyverse)
source("/project2/gilad/jpopp/sc-dynamic-eqtl/code/helpers.R")
cell.line.pca <- function(x, npc=5) {
  if (!is_tibble(x)) {
    x <- as_tibble(x)
  }
  # reshape so columns are (gene, day, and) individuals
  # filter out the days that weren't fully observed
  x <- x %>% gather(!gene, key="sample", value="counts") %>% 
            mutate(day=str_extract(sample, "[^_]+$")) %>%
            mutate(ind=str_extract(sample, "[^_]+")) %>%
            select(!sample) %>% spread(ind, counts) %>%
            drop_na()
  x <- x %>% select(!c(gene, day)) %>% apply(1, center.scale) %>% 
            t %>% as_tibble %>% 
            mutate(gene_day=paste(x$gene, x$day, sep="_")) %>%
            relocate(gene_day)
  novargenes <- x %>% select(!gene_day) %>% apply(1, var) %>% as_tibble %>%
              mutate(gene_day=x$gene_day) %>% filter(value==0) %>% .$gene_day
  x <- x %>% filter(!gene_day %in% novargenes)
  x.svd <- x %>% select(!gene_day) %>% t %>% svd(nu=npc, nv=npc)
  pcgenes <- x.svd %>% .$v %>% as_tibble %>% mutate(gene_day=x$gene_day)
  cell.line.pcs <- x.svd %>% .$u %>% as_tibble %>% 
                    `colnames<-`(paste0("PC_", seq(1, npc))) %>%
                    mutate(ind=colnames(x)[-c(1)]) %>% relocate(ind)
  var.exp <- x.svd %>% .$d %>% `^`(2) %>% `/`(sum((x.svd$d)^2))
  list("cell.line.pcs"=cell.line.pcs, "pc.genes"=pcgenes, "pve"=var.exp)
}

regular.pca <- function(x, npc=10) {
  if (!is_tibble(x)) {
    x <- as_tibble(x)
  }
  x <- x %>% select(!c(gene)) %>% apply(1, center.scale) %>% 
    t %>% as_tibble %>% mutate(gene=x$gene) %>% relocate(gene)
  novargenes <- x %>% select(!gene) %>% apply(1, var) %>% as_tibble %>%
    mutate(gene=x$gene) %>% filter(value==0) %>% .$gene
  x <- x %>% filter(!gene %in% novargenes)
  x.svd <- x %>% select(!gene) %>% t %>% svd(nu=npc, nv=npc)
  x.svd$u <- x.svd$u %>% as_tibble() %>% `colnames<-`(paste0("PC", seq(1,npc))) %>%
              mutate(sample=colnames(x)[-c(1)]) %>% relocate(sample)
  x.svd$v <- x.svd$v %>% as_tibble() %>% `colnames<-`(paste0("PC", seq(1,npc))) %>%
              mutate(gene=x$gene) %>% relocate(gene)
  x.svd
}