library(tidyverse)
library(limma)
library(pryr)
source("/project2/gilad/jpopp/sc-dynamic-eqtl/code/cell_line_pca.R")

inputArgs <- commandArgs(trailingOnly=T)

dataset <- as.character(inputArgs[1])
cis.dist <- as.character(inputArgs[2])
n.samp.pcs <- as.numeric(inputArgs[3])
n.cl.pcs <- as.numeric(inputArgs[4])
cell.type <- as.character(inputArgs[5])

message(paste0("dataset, ", dataset, 
               "\ncis.dist, ", cis.dist, 
               "\nsample pcs, ", n.samp.pcs, 
               "\ncell line pcs, ", n.cl.pcs,
               "\ncell.type, ", cell.type))

stopifnot(dataset %in% c("bulk", "pseudobulk"))
stopifnot(cis.dist %in% c("10k", "25k", "50k", "100k"))
stopifnot(n.samp.pcs<=50)
stopifnot(n.cl.pcs<=10)

# load expression, covariate data, gene-snp combinations
if (dataset == "bulk") {
  expr <- read_tsv("../data/bulk_logtpm.full.tsv")
} else if (dataset == "pseudobulk") {
  expr <- read_tsv("../data/pseudobulk_logcpm.col.full.tsv")
}
all.tests <- read_tsv(paste0("../data/gv_pairs.filtered.", cis.dist, ".tsv"))
genotypes <- read_tsv("../data/genotypes.filtered.tsv") %>%
  filter(snp %in% all.tests$snp)
colnames(genotypes) <- str_replace(colnames(genotypes), "NA", "")
genotypes <- genotypes %>% column_to_rownames("snp") %>% t %>% 
  as_tibble(rownames="ind")

# filter expression to genes we'll test and transpose 
expr <- expr %>% filter(gene %in% unique(all.tests$gene)) %>%
  column_to_rownames("gene") %>% t %>% 
  as_tibble(rownames="sample") %>%
  arrange(sample)
all.tests <- all.tests %>% filter(gene %in% colnames(expr))

# filter genotypes to snps we'll test, and expand to match expression 
inds <- expr %>% select(sample) %>% mutate(ind=str_sub(sample, 1, 5))
genotypes <- genotypes %>% select(c(ind, all.tests$snp)) %>% 
  right_join(inds, by="ind") %>% relocate(sample) %>% arrange(sample)

# get a snp list
snps <- colnames(genotypes)[-c(1,2)]

# load covariates, expand to match expression
if (dataset == "bulk") {
  cvrts <- read_tsv("../results/cibersort/bulk.inferred.tsv") %>%
    select(!c(`P-value`, Correlation, RMSE)) %>% rename(sample=Mixture) %>%
    arrange(sample) %>% select(sample, !!cell.type) 
  if (n.samp.pcs > 0) {
    samp.pcs <- read_tsv("../data/dynamic/covariates/bulk.full.pcs.tsv") %>% 
      select(c(sample, paste0("PC", seq(1, n.samp.pcs)))) %>%
      rename_with(function(s){str_replace(s, "PC", "samp_PC_")})
    cvrts <- left_join(cvrts, samp.pcs, by="sample")
  }
  if (n.cl.pcs > 0) {
    cl.pcs <- read_tsv("../data/dynamic/covariates/bulk.full.clpcs.tsv") %>% 
      mutate(ind=as.character(ind)) %>%
      right_join(inds, by="ind") %>%
      select(c(sample, paste0("PC_", seq(1, n.cl.pcs)))) %>%
      rename_with(function(s){str_replace(s, "PC_", "cl_PC_")})
    cvrts <- left_join(cvrts, cl.pcs, by="sample")
  }
} else if (dataset == "pseudobulk"){
  cvrts <- read_tsv("../results/cibersort/pseudobulk_fracs.indcol.full.true.tsv") %>% 
    arrange(sample) %>% select(sample, !!cell.type) %>% 
    filter(sample %in% expr$sample)
  if (n.samp.pcs > 0) {
    samp.pcs <- read_tsv("../data/dynamic/covariates/pseudobulk.col.full.pcs.tsv") %>% 
      select(c(sample, paste0("PC", seq(1, n.samp.pcs)))) %>%
      rename_with(function(s){str_replace(s, "PC", "samp_PC_")})
    cvrts <- left_join(cvrts, samp.pcs, by="sample")
  }
  if (n.cl.pcs > 0) {
    cl.pcs <- read_tsv("../data/dynamic/covariates/pseudobulk.day.full.clpcs.tsv") %>% 
      mutate(ind=as.character(ind)) %>%
      right_join(inds, by="ind") %>%
      select(c(sample, paste0("PC_", seq(1, n.cl.pcs)))) %>%
      rename_with(function(s){str_replace(s, "PC_", "cl_PC_")})
    cvrts <- left_join(cvrts, cl.pcs, by="sample")
  }
}

stopifnot(all.equal(expr$sample, cvrts$sample))
stopifnot(all.equal(expr$sample, genotypes$sample))

# for model specification, need cell line PC and sample PC terms in a convenient format
clpc.str <- cvrts %>%
  select(starts_with("cl_PC")) %>%
  colnames() %>%
  paste(collapse="+")
samppc.str <- cvrts %>%
  select(starts_with("samp_PC")) %>%
  colnames() %>%
  paste(collapse="+")

eqtl.call <- function(i, fname.coef, fname.res, fname.mc) {
  snp.genes = all.tests$gene[all.tests$snp==snps[i]]
  exp.snp = t(select(expr, all_of(snp.genes)))
  design = left_join(genotypes[,c(1,2,i+2)], cvrts, by="sample") %>% 
    rename(g=starts_with("chr"))
  design.snp = model.matrix(formula(paste0("~1+", samppc.str, "+(", clpc.str, "+g)*" ,cell.type)), design)
  model = lmFit(exp.snp, design.snp)
  if (anyNA(coefficients(model)[,paste0("g:", cell.type)])) {
    write_tsv(tibble(gene=snp.genes, snp=snps[i]), fname.mc, append=T)
    return()
  }
  geno.coefs = as_tibble(model$coefficients) %>% select(g)
  ct.coefs = as_tibble(model$coefficients) %>% select(cell.type)
  intx.coefs = as_tibble(model$coefficients) %>% select(paste0("g:", cell.type)) 
  intx.stdev = as_tibble(model$stdev.unscaled) %>% select(paste0("g:", cell.type))
  coef.out <- bind_cols(tibble(gene=snp.genes, snp=snps[i]), intx.coefs, intx.stdev,
                        geno.coefs, ct.coefs, tibble(sigma=as.numeric(model$sigma), df=model$df.residual))
  write_tsv(coef.out, fname.coef, append=T)
  res.out <- as_tibble(residuals(model, exp.snp, design.snp)) %>% 
    mutate(gene=snp.genes, snp=snps[i]) %>%
    relocate(gene, snp)
  write_tsv(res.out, fname.res, append=T)
}

coef.file <- paste0("../results/eqtl_dynamic/ieQTL/", dataset, "/", cell.type, "/", cis.dist, "-", n.cl.pcs, "clpcs-", n.samp.pcs, "pcs.tsv")
res.file <- paste0("../results/eqtl_dynamic/ieQTL/", dataset, "/", cell.type, "/", cis.dist, "-", n.cl.pcs, "clpcs-", n.samp.pcs, "pcs.resid.tsv")
mc.file <- paste0("../results/eqtl_dynamic/ieQTL/", dataset, "/", cell.type, "/", cis.dist, "-", n.cl.pcs, "clpcs-", n.samp.pcs, "pcs.multicollinear_tests.tsv")

if (file.exists(coef.file)) {
  file.remove(coef.file)
}
if (file.exists(res.file)) {
  file.remove(res.file)
}
if (file.exists(mc.file)) {
  file.remove(mc.file)
}
coef.header <- c("gene", "snp", paste0("g:", cell.type),
                 paste0("stdev.unscaled.g:", cell.type), 
                 "g", cell.type,
                 "sigma", "df.residual") %>% 
  purrr::map_dfc(setNames, object=list(logical())) 
write_tsv(coef.header, coef.file, col_names=T)
res.header <- c("gene", "snp", expr$sample) %>% purrr::map_dfc(setNames, object=list(logical())) 
write_tsv(res.header, res.file, col_names=T)
mc.header <- tibble("gene"=character(), "snp"=character())
write_tsv(mc.header, mc.file, col_names=T)
lapply(1:length(snps), eqtl.call, coef.file, res.file, fname.mc)
