library(tidyverse)
library(limma)
library(pryr)
source("/project2/gilad/jpopp/sc-dynamic-eqtl/code/cell_line_pca.R")
set.seed(1)

inputArgs <- commandArgs(trailingOnly=T)

dataset <- as.character(inputArgs[1])
cis.dist <- as.character(inputArgs[2])
n.samp.pcs <- as.numeric(inputArgs[3])
n.cl.pcs <- as.numeric(inputArgs[4])
agg <- as.character(inputArgs[5])

message(paste0("dataset, ", dataset, 
               "\ncis.dist, ", cis.dist, 
               "\nsample pcs, ", n.samp.pcs,
               "\ncell line pcs, ", n.cl.pcs,
               "\naggregate, ", agg))

stopifnot(dataset %in% c("bulk", "bulk7", "pseudobulk", "pseudobulk-cm", "pseudobulk-cf"))
stopifnot(cis.dist %in% c("10k", "25k", "50k", "100k"))
stopifnot(n.samp.pcs<=50)
stopifnot(n.cl.pcs<=10)

# load expression, covariate data, gene-snp combinations
if (dataset == "bulk") {
  expr <- read_tsv("../data/bulk_logtpm.full.tsv")
} else if (dataset == "bulk7") {
  days.inc <- c("_0", "_1", "_3", "_5", "_7", "_11", "_15")
  expr <- read_tsv("../data/bulk_logtpm.full.tsv") %>% 
    select(c(gene, ends_with(days.inc)))
} else if (dataset == "pseudobulk") {
  expr <- read_tsv(paste0("../data/pseudobulk_logcpm.", agg, ".full.tsv"))
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

# get time points
if ("day" %in% agg) {
  days <- expr %>% select(sample) %>% mutate(day=as.numeric(str_sub(sample, 7)))
} else {
  medians <- read_tsv(paste0("../results/eqtl_dynamic/linear_dQTL/pseudobulk/", agg, "/15bin_medians.tsv")) %>%
    rename(sample=bin.ind)
  days <- expr %>% select(sample) %>% left_join(medians, by="sample") %>%
    rename(day=t)
}

# get a snp list
snps <- colnames(genotypes)[-c(1,2)] %>%
  as_tibble() %>%
  rename(snp=value) %>%
  rowid_to_column("snp_num")
all.tests <- right_join(all.tests, snps, by="snp")


# load covariates, expand to match expression
if (dataset == "bulk") {
  cvrts <- read_tsv("../data/dynamic/covariates/bulk.full.clpcs.tsv") %>% 
    select(1:(n.cl.pcs+1)) %>% mutate(ind=as.character(ind)) %>% 
    right_join(inds, by="ind") %>% select(!ind) %>%
    arrange(sample) %>%
    rename_with(function(s){str_replace(s, "PC_", "cl_PC_")})
  if (n.samp.pcs > 0) {
    samp.pcs <- read_tsv("../data/dynamic/covariates/bulk.full.pcs.tsv") %>%
      select(1:(n.samp.pcs+1)) %>%
      rename_with(function(s){str_replace(s,"PC","samp_PC_")})
    cvrts <- left_join(cvrts, samp.pcs, by="sample")
  }
} else if (dataset == "bulk7") {
  cvrts <- read_tsv("../data/dynamic/covariates/bulk.full.clpcs.tsv") %>% 
    select(1:(n.cl.pcs+1)) %>% mutate(ind=as.character(ind)) %>% 
    right_join(inds, by="ind") %>% select(!ind) %>%
    arrange(sample) %>%
    rename_with(function(s){str_replace(s, "PC_", "cl_PC_")})
  if (n.samp.pcs > 0) {
    samp.pcs <- read_tsv("../data/dynamic/covariates/bulk.full.pcs.tsv") %>%
      select(1:(n.samp.pcs+1)) %>%
      rename_with(function(s){str_replace(s,"PC","samp_PC_")})
    cvrts <- left_join(cvrts, samp.pcs, by="sample")
  }
} else if (dataset == "pseudobulk"){
  cvrts <- read_tsv(paste0("../data/dynamic/covariates/pseudobulk.", agg, ".full.clpcs.tsv")) %>% 
    select(1:(n.cl.pcs+1))%>% mutate(ind=as.character(ind)) %>% 
    right_join(inds, by="ind") %>% select(!ind) %>%
    arrange(sample) %>%
    rename_with(function(s){str_replace(s, "PC_", "cl_PC_")})
  if (n.samp.pcs > 0) {
    samp.pcs <- read_tsv(paste0("../data/dynamic/covariates/pseudobulk.", agg, ".full.pcs.tsv")) %>%
      select(1:(n.samp.pcs+1)) %>%
      rename_with(function(s){str_replace(s,"PC","samp_PC_")})
    cvrts <- left_join(cvrts, samp.pcs, by="sample")
  }
}

stopifnot(sum(expr$sample != genotypes$sample) == 0)
stopifnot(sum(expr$sample != days$sample) == 0)

# for model specification, need cell line PC and sample PC terms in a convenient format
clpc.str <- cvrts %>%
  select(starts_with("cl_PC")) %>%
  colnames() %>%
  paste(collapse="+")
samppc.str <- cvrts %>%
  select(starts_with("samp_PC")) %>%
  colnames() %>%
  paste(collapse="+")
if (n.samp.pcs > 0 & n.cl.pcs > 0) {
  formula.str <- paste0("~1+", samppc.str, "+(", clpc.str, ")*day+g+g:day.permute")
} else if (n.samp.pcs > 0) {
  formula.str <- paste0("~1+", samppc.str, "+day+g+g:day.permute")
} else if (n.cl.pcs > 0) {
  formula.str <- paste0("~1+(", clpc.str, ")*day+g+g:day.permute")
} else {
  formula.str <- paste0("~1+g+day+g:day.permute")
}

eqtl.call <- function(i, fname.coef) {
  genei = all.tests$gene[i]
  snpi = all.tests$snp[i]
  expi = t(select(expr, genei))
  design = left_join(genotypes[,c(1,2,all.tests$snp_num[i]+2)], days, by="sample") %>% 
    rename(g=starts_with("chr")) %>%
    left_join(cvrts, by="sample") %>%
    mutate(day.permute=sample(day, nrow(.)))
  design.snp = model.matrix(formula(formula.str), design)
  model=lmFit(expi, design.snp)
  if (anyNA(coefficients(model)[,"g:day.permute"])) {
    write_tsv(tibble(gene=genei, snp=snpi), fname.mc, append=T)
    return()
  }
  coef.out = tibble("gene"=genei, "snp"=snpi, "beta.g"=as.numeric(model$coefficients[,"g"]),
                    "beta.gxt"=as.numeric(model$coefficients[,"g:day.permute"]),
                    "beta.t"=as.numeric(model$coefficients[,"day"]),
                    "stdev.unscaled.g"=as.numeric(model$stdev.unscaled[,"g"]), 
                    "stdev.unscaled.gxt"=as.numeric(model$stdev.unscaled[,"g:day.permute"]),
                    "stdev.unscaled.t"=as.numeric(model$stdev.unscaled[,"day"]),
                    "sigma"=as.numeric(model$sigma), "df.residual"=as.numeric(model$df.residual))
  write_tsv(coef.out, fname.coef, append=T)
}

coef.file <- paste0("../results/eqtl_dynamic/linear_dQTL/", dataset, "/", agg, "/", cis.dist, "-", n.cl.pcs, "clpcs-", n.samp.pcs, "pcs.permuted.tsv")

if (file.exists(coef.file)) {
  file.remove(coef.file)
}
  
coef.header <- c("gene", "snp", "beta.g", "beta.gxt", "beta.t",
                 "stdev.unscaled.g", "stdev.unscaled.gxt", "stdev.unscaled.t",
                  "sigma", "df.residual") %>% 
  purrr::map_dfc(setNames, object=list(logical())) 
write_tsv(coef.header, coef.file, col_names=T)

lapply(1:nrow(all.tests), eqtl.call, coef.file)
