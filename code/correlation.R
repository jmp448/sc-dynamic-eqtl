library(limma)
library(tidyverse)

## regress out PCs for single time points
for (dataset in c("bulk", "pseudobulk")) {
  for (npcs in c(1, 2, 3, 4, 5)) {
    for (d in paste0("day", c(0, 1, 3, 5, 7, 11, 15))) {
      # load expression and pcs
      pcs <- read_tsv(paste0("../data/static/", dataset, "/day/", d, "/covariates.tsv")) %>% select(1:(1+npcs))
      exp <- read_tsv(paste0("../data/static/", dataset, "/day/", d, "/expression.tsv")) %>%
        column_to_rownames("sample") %>% t %>% as_tibble(rownames="gene")
      stopifnot(all.equal(colnames(exp)[-c(1)], as.character(pcs$sample)))
      
      # use lmFit to quickly regress out the specified number of expression PCs
      design <- model.matrix(~.-sample, pcs)
      response <- exp %>% select(!gene)
      model <- lmFit(response, design)
      
      # save the residuals
      resids <- resid(model, response) %>% as_tibble %>% mutate(gene=exp$gene) %>%
        column_to_rownames("gene") %>% t %>% as_tibble(rownames="sample") %>%
        write_tsv(paste0("../data/static/", dataset, "/day/", d, "/expression_reg", npcs, ".tsv"))
    }
  }
}

## regress out cell line PCs for all time points
for (dataset in c("bulk", "pseudobulk")) {
  for (npcs in c(1, 2, 3, 4, 5)) {
    # load expression and pcs
    pcs <- read_tsv(paste0("../data/static/", dataset, "/day/", d, "/covariates.tsv")) %>% select(1:(1+npcs))
    exp <- read_tsv(paste0("../data/static/", dataset, "/day/", d, "/expression.tsv")) %>%
      column_to_rownames("sample") %>% t %>% as_tibble(rownames="gene")
    stopifnot(all.equal(colnames(exp)[-c(1)], as.character(pcs$sample)))
    
    # use lmFit to quickly regress out the specified number of expression PCs
    design <- model.matrix(~.-sample, pcs)
    response <- exp %>% select(!gene)
    model <- lmFit(response, design)
    
    # save the residuals
    resids <- resid(model, response) %>% as_tibble %>% mutate(gene=exp$gene) %>%
      column_to_rownames("gene") %>% t %>% as_tibble(rownames="sample") %>%
      write_tsv(paste0("../data/static/", dataset, "/day/", d, "/expression_reg", npcs, ".tsv"))
  }
}

# measure correlation across all genes between individuals
cors <- tibble("rho"=numeric(), "cl1"=character(), "cl2"=character(), "npcs"=numeric(), "day"=character())
for (npc in seq(0, 5)) {
  for (d in paste0("day", c(0, 1, 3, 5, 7, 11, 15))) {
    if (npc == 0) {
      bulk <- read_tsv(paste0("../data/static/bulk/day/", d, "/expression.tsv")) 
      pseudobulk <- read_tsv(paste0("../data/static/pseudobulk/day/", d, "/expression.tsv")) 
    } else if (npc > 0) {
      bulk <- read_tsv(paste0("../data/static/bulk/day/", d, "/expression_reg", npc, ".tsv")) 
      pseudobulk <- read_tsv(paste0("../data/static/pseudobulk/day/", d, "/expression_reg", npc, ".tsv")) 
    }
    bulk <- bulk %>% 
      column_to_rownames("sample") %>%
      t %>% as_tibble(rownames="gene")
    pseudobulk <- pseudobulk %>% 
      column_to_rownames("sample") %>%
      t %>% as_tibble(rownames="gene")
    gene.set <- intersect(bulk$gene, pseudobulk$gene)
    bulk <- bulk %>% filter(gene %in% gene.set) %>% arrange(gene)
    pseudobulk <- pseudobulk %>% filter(gene %in% gene.set) %>% arrange(gene)
    inds1 <- c()
    inds2 <- c()
    for (i in colnames(bulk)[-c(1)]) {
      for (j in colnames(pseudobulk)[-c(1)]) {
        inds1 <- c(inds1, i)
        inds2 <- c(inds2, j)
      }
    }
    cors.ij <- mapply(function(i1, i2){cor(bulk[[i1]], pseudobulk[[i2]])}, inds1, inds2) 
    cors <- bind_rows(cors, tibble("rho"=cors.ij, "cl1"=inds1, "cl2"=inds2, "npcs"=npc, "day"=d))
  }
}
write_tsv(cors, "../results/eqtl_static/correlation.tsv")

# measure correlation across all individuals for each gene
# null will be permutations
cor.comp <- function(d, npc) {
  if (npc == 0) {
    bulk <- read_tsv(paste0("../data/static/bulk/day/", d, "/expression.tsv")) 
    pseudobulk <- read_tsv(paste0("../data/static/pseudobulk/day/", d, "/expression.tsv"))
  } else if (npc > 0) {
    bulk <- read_tsv(paste0("../data/static/bulk/day/", d, "/expression_reg", npc, ".tsv")) 
    pseudobulk <- read_tsv(paste0("../data/static/pseudobulk/day/", d, "/expression_reg", npc, ".tsv")) 
  }
  bulk <- bulk %>% select(intersect(colnames(bulk), colnames(pseudobulk)))
  pseudobulk <- pseudobulk %>% select(intersect(colnames(bulk), colnames(pseudobulk)))
  null1 <- pseudobulk %>% sample_frac(1)
  null2 <- pseudobulk %>% sample_frac(1)
  cors.same <- tibble("rho"=sapply(colnames(bulk)[-c(1)], function(g){cor(bulk[[g]], pseudobulk[[g]])}),
                      "permuted"=F,
                      "day"=d,
                      "npcs"=npc)
  cors.diff <- tibble("rho"=c(
    sapply(colnames(bulk)[-c(1)], function(g){cor(bulk[[g]], null1[[g]])}),
    sapply(colnames(bulk)[-c(1)], function(g){cor(bulk[[g]], null2[[g]])})),
    "permuted"=T,
    "day"=d,
    "npcs"=npc)
  cors <- bind_rows(cors.same, cors.diff)
  cors
}
pcs <- c()
days <- c()
for (npc in seq(0, 3)) {
  for (d in paste0("day", c(0, 1, 3, 5, 7, 11, 15))) {
    pcs <- c(pcs, npc)
    days <- c(days, d)
  }
}
all.cors <- bind_rows(pmap_dfr(list(days, pcs), cor.comp))
write_tsv(all.cors, "../results/eqtl_static/correlation2.tsv")