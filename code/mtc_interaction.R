library(tidyverse)
library(stats)
library(qvalue)
library(mashr)
library(limma)

inputArgs <- commandArgs(trailingOnly=T)

dataset <- as.character(inputArgs[1])
cis.dist <- as.character(inputArgs[2])
n.samp.pcs <- as.numeric(inputArgs[3])
n.cl.pcs <- as.numeric(inputArgs[4])
cell.type <- as.character(inputArgs[5])

message(paste0("dataset, ", dataset, 
               "\ncis.dist, ", cis.dist, 
               "\nsamp pcs, ", n.samp.pcs, 
               "\ncl pcs, ", n.cl.pcs,
               "\ncell.type, ", cell.type))

stopifnot(dataset %in% c("bulk", "pseudobulk"))
stopifnot(cis.dist %in% c("10k", "25k", "50k", "100k"))
stopifnot(n.samp.pcs<=50)
stopifnot(n.cl.pcs<=20)

results <- read_tsv(paste0("../results/eqtl_dynamic/ieQTL/", dataset, "/", cell.type, "/", cis.dist, "-", n.cl.pcs, "clpcs-", n.samp.pcs ,"pcs.tsv"))

results.ct <- results %>%
  rename(stdev.unscaled=paste0("stdev.unscaled.g:", cell.type)) %>%
  rename(coefficients=paste0("g:", cell.type)) %>%
  mutate(se.unadj = stdev.unscaled*sigma) %>%
  mutate(t.unadj = coefficients / se.unadj) %>%
  mutate(p.unadj = 2*pt(-abs(t.unadj), df=df.residual[1]))

# Bonferroni FWER Correction
gene_ct <- table(results.ct$gene)
results.ct$bonf.p.unadj <- sapply(gene_ct[results.ct$gene] * results.ct$p.unadj, function(x){min(1,x)})

# perform EB shrinkage on the standard errors
results.ct <- results.ct %>% eBayes %>%
  rename(t.adj=t) %>%
  mutate(p.adj=map_dbl(t.adj, function(t, nu=df.total[1])(2*pt(-abs(t), df=nu))))
results.ct$bonf.p.adj <- sapply(gene_ct[results.ct$gene] * results.ct$p.adj, function(x){min(1,x)})

# run ashr
betahat <- results.ct$coefficients
sehat <- results.ct$stdev.unscaled * sqrt(results.ct$s2.post)
df <- results.ct$df.total[1]
ashr.fit <- ash(betahat, sehat, df=df)

results.ct <- results.ct %>%
  mutate(lfsr=ashr.fit$result$lfsr)

# save the fitted model
g <- get_fitted_g(ashr.fit)
saveRDS(g, paste0("../results/eqtl_dynamic/ieQTL/", dataset, "/", cell.type, "/", cis.dist, "-", n.cl.pcs, "clpcs-", n.samp.pcs ,"pcs.ashr.g.rds"))

# subset to top hits per gene
top.hits  <- results.ct %>% arrange(p.adj) %>%
  filter(!duplicated(gene)) %>%
  mutate(qval.unadj=qvalue(bonf.p.unadj)$q) %>%
  mutate(qval.adj=qvalue(bonf.p.adj)$q)

write_tsv(results.ct, paste0("../results/eqtl_dynamic/ieQTL/", dataset, "/", cell.type, "/", cis.dist, "-", n.cl.pcs, "clpcs-", n.samp.pcs ,"pcs.adj.tsv"))
write_tsv(top.hits, paste0("../results/eqtl_dynamic/ieQTL/", dataset, "/", cell.type, "/", cis.dist, "-", n.cl.pcs, "clpcs-", n.samp.pcs ,"pcs.tophits.tsv"))
