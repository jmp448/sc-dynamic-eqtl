library(tidyverse)
library(limma)

inputArgs <- commandArgs(trailingOnly=T)

dataset <- as.character(inputArgs[1])
aggregate <- as.character(inputArgs[2])
cis.dist <- as.character(inputArgs[3])
npcs <- as.numeric(inputArgs[4])
cell.type <- as.character(inputArgs[5])

expr <- read_tsv(paste0("../data/static/", dataset, "/", aggregate, "/", cell.type, "/expression.tsv")) %>%
  mutate(ind=as.character(sample), .keep="unused", .after=1)
covariates <- read_tsv(paste0("../data/static/", dataset, "/", aggregate, "/", cell.type, "/covariates.tsv")) %>% 
  select(1:(1+npcs)) %>%
  mutate(ind=as.character(sample), .keep="unused", .after=1)
genotypes <- read_tsv("../data/genotypes.filtered.tsv")
colnames(genotypes) <- colnames(genotypes) %>% str_replace("NA", "")
genotypes <- genotypes %>% select(c(snp, expr$ind))

all.tests <- read_tsv(paste0("../data/gv_pairs.filtered.", cis.dist, ".tsv")) %>%
  filter(gene %in% colnames(expr))

# filter expression to genes we'll test
expr <- select(expr, c(ind, unique(all.tests$gene)))

# filter genotypes to snps we'll test
genotypes <- filter(genotypes, snp %in% all.tests$snp)
snps <- genotypes$snp
genotypes <- select(genotypes, !snp)
stopifnot(nrow(genotypes)==length(snps))
# reorder genotypes to match expr
stopifnot(nrow(expr) == ncol(genotypes))
genotypes <- relocate(genotypes, expr$ind)
stopifnot(sum(colnames(genotypes)!=expr$ind)==0)

# make sure expr and covariate order is the same
stopifnot(all.equal(expr$ind, covariates$ind))

eqtl.call <- function(i) {
  snp.genes = all.tests$gene[all.tests$snp==snps[i]]
  exp = t(select(expr, all_of(snp.genes)))
  design = mutate(covariates, genotype=as.numeric(genotypes[i,]), intercept=1)
  design = select(relocate(design, intercept), !ind)
  model = lmFit(exp, design)
  tibble("gene"=snp.genes, "snp"=snps[i], "coefficients"=as.numeric(model$coefficients[,"genotype"]),
                 "stdev.unscaled"=as.numeric(model$stdev.unscaled[,"genotype"]), 
                 "sigma"=as.numeric(model$sigma), "df.residual"=as.numeric(model$df.residual))
}

results <- bind_rows(lapply(1:length(snps), eqtl.call))
write_tsv(results, paste0("../results/eqtl_static/", dataset, "/", aggregate, "/", cell.type, "/", cis.dist, "-", npcs, "pcs.tsv"))