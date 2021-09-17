library(tidyverse)
library(limma)
library(vroom)

inputArgs <- commandArgs(trailingOnly=T)
# called by ../wrappers/eqtl_static.sh, params from ../wrappers/eqtl_static_queue.sh

dataset <- as.character(inputArgs[1])
aggregate <- as.character(inputArgs[2])
cis.dist <- as.character(inputArgs[3])
npcs <- as.numeric(inputArgs[4])
cell.type <- as.character(inputArgs[5])
chunk_i <- as.numeric(inputArgs[6])
chunks_tot <- as.numeric(inputArgs[7])

expr <- read_tsv(paste0("../data/static/", dataset, "/", aggregate, "/", cell.type, "/expression.tsv")) %>%
  mutate(ind=as.character(sample), .keep="unused", .after=1)
covariates <- read_tsv(paste0("../data/static/", dataset, "/", aggregate, "/", cell.type, "/pcs.tsv")) %>% 
  select(1:(1+npcs)) %>%
  mutate(ind=as.character(sample), .keep="unused", .after=1)
all.tests <- read_tsv(paste0("../data/static/", dataset, "/", aggregate, "/", cell.type, "/filtered_tests.", cis.dist, ".tsv")) %>%
  arrange(snp)

# subset to a fraction (1/chunks_tot) of the total snps
snps <- unique(all.tests$snp)
snps_i <- split(snps, cut(seq_along(snps), chunks_tot, labels = FALSE))[[chunk_i]]
all.tests <- filter(all.tests, snp %in% snps_i)

# filter expression to genes we'll test
expr <- select(expr, c(ind, unique(all.tests$gene)))

# filter genotypes to snps we'll test
genotypes <- vroom("../data/genotypes.tsv") %>%
  filter(snp %in% all.tests$snp) %>% 
  select(c(snp, expr$ind))
snps <- genotypes$snp
genotypes <- select(genotypes, !snp)
stopifnot(nrow(genotypes)==length(snps))

# reorder genotypes to match expr
stopifnot(nrow(expr) == ncol(genotypes))
genotypes <- relocate(genotypes, expr$ind)
stopifnot(sum(colnames(genotypes)!=expr$ind)==0)

# make sure expr and covariate order is the same
stopifnot(all.equal(expr$ind, covariates$ind))

eqtl.call <- function(i, fout) {
  snp.genes = all.tests$gene[all.tests$snp==snps[i]]
  exp = t(select(expr, all_of(snp.genes)))
  design = mutate(covariates, genotype=as.numeric(genotypes[i,]), intercept=1)
  design = select(relocate(design, intercept), !ind)
  model = lmFit(exp, design)
  write_tsv(tibble("gene"=snp.genes, "snp"=snps[i], "coefficients"=as.numeric(model$coefficients[,"genotype"]),
                 "stdev.unscaled"=as.numeric(model$stdev.unscaled[,"genotype"]), 
                 "sigma"=as.numeric(model$sigma), "df.residual"=as.numeric(model$df.residual)),
            fout, append=T)
}

# if the output file exists, overwrite it
out.file <- paste0("../results/eqtl_static/", dataset, "/", aggregate, "/", cell.type, "/", cis.dist, "-", npcs, "pcs-", chunk_i, ".tsv")
if (file.exists(out.file)) {
  file.remove(out.file)
}

# create header (outside of loop)
out.header <- c("gene", "snp", "coefficients", "stdev.unscaled", "sigma", "df.residual") %>% 
  purrr::map_dfc(setNames, object=list(logical())) %>%
  write_tsv(out.file, col_names=T)

system.time(lapply(1:length(snps), eqtl.call, out.file))
