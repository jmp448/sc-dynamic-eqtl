library(tidyverse)
library(vroom)
library(limma)
library(pryr)
source("/project2/gilad/jpopp/sc-dynamic-eqtl/code/cell_line_pca.R")
set.seed(42)

inputArgs <- commandArgs(trailingOnly=T)

dataset <- as.character(inputArgs[1])
cis.dist <- as.character(inputArgs[2])
n.samp.pcs <- as.numeric(inputArgs[3])
n.cl.pcs <- as.numeric(inputArgs[4])
cell.type.reg <- as.logical(inputArgs[5])
cell.type <- as.character(inputArgs[6])
chunk_i <- as.numeric(inputArgs[7])
chunks_tot <- as.numeric(inputArgs[8])

message(paste0("dataset, ", dataset, 
               "\ncis.dist, ", cis.dist, 
               "\nsample pcs, ", n.samp.pcs, 
               "\ncell line pcs, ", n.cl.pcs,
               "\ncell type regression, ", cell.type.reg,
               "\ncell.type, ", cell.type,
               "\nchunk index, ", chunk_i,
               "\nchunks_tot, ", chunks_tot))

stopifnot(dataset %in% c("bulk", "bulk7", "pseudobulk"))
stopifnot(cis.dist %in% c("10k", "25k", "50k", "100k"))
stopifnot(n.samp.pcs<=50)
stopifnot(n.cl.pcs<=10)

# load expression, covariate data, gene-snp combinations
if (dataset %in% c("bulk", "bulk7")) {
  expr <- read_tsv(paste0("../data/", dataset, "/day/logtpm.tsv")) 
} else {
  expr <- read_tsv(paste0("../data/", dataset, "/day/logcpm.tsv")) 
}
expr <- expr %>% column_to_rownames("gene") %>% apply(1, center.scale) %>% t %>% as_tibble(rownames="gene") # standardized expression
all.tests <- read_tsv(paste0("../data/", dataset, "/day/filtered_tests.", cis.dist, ".tsv")) %>%
  arrange(snp)

# subset to a fraction (1/chunks_tot) of the total snps
snps <- unique(all.tests$snp)
snps_i <- split(snps, cut(seq_along(snps), chunks_tot, labels = FALSE))[[chunk_i]]
all.tests <- filter(all.tests, snp %in% snps_i)

# filter expression to genes we'll test and transpose 
expr <- expr %>% filter(gene %in% unique(all.tests$gene)) %>%
  column_to_rownames("gene") %>% t %>% 
  as_tibble(rownames="sample") %>%
  arrange(sample)
all.tests <- all.tests %>% 
  filter(gene %in% colnames(expr))

# filter genotypes to snps we'll test
genotypes <- vroom("../data/genotypes.tsv") %>%
  filter(snp %in% all.tests$snp)

# impose the sample repeat structure into the genotype matrix
genotypes <- genotypes %>% column_to_rownames("snp") %>% t %>% 
  as_tibble(rownames="ind") %>% arrange(ind)
inds <- expr %>% dplyr::select(sample) %>% mutate(ind=str_sub(sample, 1, 5)) %>%
  arrange(ind)
ind.indices <- as.numeric(factor(inds$ind))
geno.mat <- genotypes %>% column_to_rownames("ind") %>% as.matrix
genotypes <- geno.mat[ind.indices,] %>% `rownames<-`(inds$sample) %>%
  as_tibble(rownames="sample")

# get a new snp list
snps <- colnames(genotypes)[-c(1)]

# load covariates, expand to match expression
cvrts <- read_tsv(paste0("../data/", dataset, "/day/clpcs.tsv")) %>%
  select(1:(n.cl.pcs+1)) %>% mutate(ind=as.character(ind)) %>% 
  right_join(inds, by="ind") %>% select(!ind) %>%
  rename_with(function(s){str_replace(s, "PC_", "cl_PC_")})
if (n.samp.pcs > 0) {
  samp.pcs <- read_tsv(paste0("../data/", dataset, "/day/pcs.tsv")) %>%
    select(1:(n.samp.pcs+1)) %>%
    rename_with(function(s){str_replace(s,"PC","samp_PC_")}) %>%
    mutate(ind=str_extract(sample, "[^_]+"))
  cvrts <- right_join(cvrts, samp.pcs, by="sample") 
}
if (dataset == "bulk") {
  ct.props <- read_tsv("../results/cibersort/bulk.inferred.tsv") %>%
    select(!c(`P-value`, Correlation, RMSE)) %>% rename(sample=Mixture) %>%
    arrange(sample) %>% select(!UNK)
} else if (dataset == "bulk7") {
  ct.props <- read_tsv("../results/cibersort/bulk.inferred.tsv") %>%
    select(!c(`P-value`, Correlation, RMSE)) %>% rename(sample=Mixture) %>%
    arrange(sample) %>% select(!UNK) %>% filter(sample %in% expr$sample)
} else {
  ct.props <- read_tsv(paste0("../results/cibersort/", dataset, "_indday.full.true.tsv")) %>%
    arrange(sample) %>% select(!UNK)
}
cell.types <- setdiff(colnames(ct.props), "sample")
cvrts <- left_join(cvrts, ct.props, by="sample")

stopifnot(all.equal(expr$sample, genotypes$sample))
stopifnot(all.equal(expr$sample, cvrts$sample))

# for model specification, need cell line PC and sample PC terms in a convenient format
selective.append <- function(s) {
  if ("" %in% s) {
    s
  } else {
    append(s, "")
  }
}
clpc.str <- cvrts %>%
  select(starts_with("cl_PC")) %>%
  colnames() %>%
  selective.append %>%
  paste(collapse="+")
samppc.str <- cvrts %>%
  select(starts_with("samp_PC")) %>%
  colnames() %>%
  selective.append %>%
  paste(collapse="+")
if (cell.type.reg) {
  ctprop.str <- setdiff(c("IPSC", "MES", "CMES", "PROG", "CM", "CF"), cell.type) %>%
    selective.append %>%
    paste(collapse="+")
} else {
  ctprop.str <- ""
}

formula.str <- paste0("~1+", samppc.str, ctprop.str, "(", clpc.str, "g)*", cell.type)

eqtl.call <- function(i, fname.coef, fname.res, fname.mc) {
  snp.genes = all.tests$gene[all.tests$snp==snps[i]]
  exp.snp = t(select(expr, all_of(snp.genes)))
  design = left_join(genotypes[,c(1,i+1)], cvrts, by="sample") %>% 
    rename(g=starts_with("chr"))
  design.snp = model.matrix(formula(formula.str), design)
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

if (cell.type.reg) {
  cell.type.reg <- "regtypes"
} else {
  cell.type.reg <- "notypes"
}
coef.file <- paste0("../results/eqtl_dynamic/ieQTL/", dataset, "/", cell.type, "/", cis.dist, "-", n.cl.pcs, "clpcs-", n.samp.pcs, "pcs-", cell.type.reg, "-", chunk_i, ".tsv")
res.file <- paste0("../results/eqtl_dynamic/ieQTL/", dataset, "/", cell.type, "/", cis.dist, "-", n.cl.pcs, "clpcs-", n.samp.pcs, "pcs-", cell.type.reg, "-", chunk_i, ".resid.tsv")
mc.file <- paste0("../results/eqtl_dynamic/ieQTL/", dataset, "/", cell.type, "/", cis.dist, "-", n.cl.pcs, "clpcs-", n.samp.pcs, "pcs-", cell.type.reg, "-", chunk_i, ".multicollinear_tests.tsv")

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
