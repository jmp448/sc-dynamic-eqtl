library(tidyverse)
library(mashr)
set.seed(1234)

inputArgs <- commandArgs(trailingOnly=T)

dataset <- as.character(inputArgs[1])
agg <- as.character(inputArgs[2])
cis.dist <- as.character(inputArgs[3])
npcs <- as.numeric(inputArgs[4])

get_coefs <- function(t) {
  results.t <- read_tsv(paste0("../results/eqtl_static/", dataset, "/", agg, "/", t, "/", cis.dist, "-", npcs, "pcs.mtc.tsv"))
  beta.hat <- results.t %>% 
    mutate("{t}":=coefficients) %>%
    mutate(gv=paste(gene, snp, sep="--")) %>%
    select(gv, !!t)
  beta.hat
} 

get_se <- function(t) {
  results.t <- read_tsv(paste0("../results/eqtl_static/", dataset, "/", agg, "/", t, "/", cis.dist, "-", npcs, "pcs.mtc.tsv"))
  se.hat <- results.t %>% 
    mutate("{t}":=se.adj) %>%
    mutate(gv=paste(gene, snp, sep="--")) %>%
    select(gv, !!t)
  se.hat
} 

get_dfs <- function(t) {
  results.t <- read_tsv(paste0("../results/eqtl_static/", dataset, "/", agg, "/", t, "/", cis.dist, "-", npcs, "pcs.mtc.tsv"))
  df.total <- results.t %>% 
    mutate("{t}":=df.total) %>%
    mutate(gv=paste(gene, snp, sep="--")) %>%
    select(gv, !!t)
  df.total
} 
 
if (dataset == "bulk") {
  days <- paste0("day", seq(0, 15))
} else if ((dataset == "pseudobulk") & (agg == "day")) {
  days <- paste0("day", c(0, 1, 3, 5, 7, 11, 15))
} else if ((dataset == "pseudobulk") & (agg == "type")) {
  days <- c("iPSC", "meso", "EMT", "cardiomes", "prog", "CM", "EPDC")
}

beta.hat <- map(days, get_coefs) %>% reduce(inner_join, by="gv") %>%
  column_to_rownames("gv") %>% as.matrix
se.hat <- map(days, get_se) %>% reduce(inner_join, by="gv") %>%
  column_to_rownames("gv") %>% as.matrix
df.adj <- map(days, get_dfs) %>% reduce(inner_join, by="gv") %>%
  column_to_rownames("gv") %>% as.matrix

# identify top hit per gene across all types/days
top.hits <- beta.hat %>% `/`(se.hat) %>% abs %>% 
  apply(1, max) %>%
  as_tibble(rownames="gv") %>%
  mutate(gene=str_extract(gv, "[^-]+")) %>%
  arrange(desc(value)) %>%
  filter(!duplicated(gene)) %>%
  .$gv

# identify a random subset of 250000 tests
random.hits = sample(rownames(beta.hat), min(250000, dim(beta.hat)[1]))

# estimate null correlation
data.temp = mash_set_data(beta.hat[random.hits,],se.hat[random.hits,])
Vhat = estimate_null_correlation_simple(data.temp)
rm(data.temp)

data.random = mash_set_data(beta.hat[random.hits,],se.hat[random.hits,],V=Vhat, df=df.adj[random.hits,])
data.strong = mash_set_data(beta.hat[top.hits,],se.hat[top.hits,],V=Vhat, df=df.adj[top.hits,])

# get data-driven covariances
U.pca = cov_pca(data.strong,5)
U.ed = cov_ed(data.strong, U.pca)

# fit mash model
U.c = cov_canonical(data.random)
m = mash(data.random, Ulist = c(U.ed,U.c), outputlevel = 1)

m2 = mash(data.strong, g=get_fitted_g(m), fixg=TRUE)

saveRDS(m2, paste0("../results/eqtl_static/", dataset, "/", agg, "/", cis.dist, "-", npcs, "pcs.mashr.tsv"))