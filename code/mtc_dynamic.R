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
agg <- as.character(inputArgs[5])

message(paste0("dataset, ", dataset, 
               "\ncis.dist, ", cis.dist, 
               "\nsample pcs, ", n.samp.pcs,
               "\ncell line pcs, ", n.cl.pcs,
               "\naggregate, ", agg))

stopifnot(dataset %in% c("bulk", "bulk7", "pseudobulk", "pseudobulk-drop2"))
stopifnot(cis.dist %in% c("10k", "25k", "50k", "100k"))
stopifnot(n.samp.pcs<=50)
stopifnot(n.cl.pcs<=10)
stopifnot(agg %in% c("day", "cmbin", "epdcbin"))

results <- read_tsv(paste0("../results/eqtl_dynamic/linear_dQTL/", 
                           dataset, "/", agg, "/", cis.dist, "-", n.cl.pcs, "clpcs-", n.samp.pcs, "pcs.tsv"))

# Bonferroni multiple testing correction
results <- results %>%
  mutate(stdev.unscaled=stdev.unscaled.gxt) %>%
  mutate(coefficients=beta.gxt) %>%
  mutate(se.unadj = stdev.unscaled.gxt*sigma) %>%
  mutate(t.unadj = beta.gxt / se.unadj) %>%
  mutate(p.unadj = 2*pt(-abs(t.unadj), df=df.residual[1]))

gene_ct <- table(results$gene)
results$bonf.p.unadj <- sapply(gene_ct[results$gene] * results$p.unadj, function(x){min(1,x)})

# perform EB shrinkage on the standard errors
results <- results %>% eBayes %>%
  rename(t.adj=t) %>%
  mutate(p.adj=map_dbl(t.adj, function(t, nu=df.total[1])(2*pt(-abs(t), df=nu))))
results$bonf.p.adj <- sapply(gene_ct[results$gene] * results$p.adj, function(x){min(1,x)})

# run ashr
betahat <- results$beta.gxt
sehat <- results$stdev.unscaled.gxt * sqrt(results$s2.post)
df <- results$df.total[1]
ashr.fit <- ash(betahat, sehat, df=df)

results <- results %>%
  mutate(lfsr=ashr.fit$result$lfsr)

# save the fitted model
g <- get_fitted_g(ashr.fit)
saveRDS(g, paste0("../results/eqtl_dynamic/linear_dQTL/", dataset, "/", agg, "/", cis.dist, "-", n.cl.pcs, "clpcs-", n.samp.pcs, "pcs.ashr.g.rds"))

# subset to top hits per gene and compute Storey's q value
top.hits  <- results %>% arrange(p.unadj) %>%
  filter(!duplicated(gene)) %>%
  mutate(qval.unadj=qvalue(bonf.p.unadj)$q) %>%
  mutate(qval.adj=qvalue(bonf.p.adj)$q)

write_tsv(results, paste0("../results/eqtl_dynamic/linear_dQTL/", dataset, "/", agg, "/", cis.dist, "-", n.cl.pcs, "clpcs-", n.samp.pcs, "pcs.mtc.tsv"))
write_tsv(top.hits, paste0("../results/eqtl_dynamic/linear_dQTL/", dataset, "/", agg, "/", cis.dist, "-", n.cl.pcs, "clpcs-", n.samp.pcs, "pcs.tophits.tsv"))
