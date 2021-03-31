library(tidyverse)
library(stats)
library(qvalue)
library(mashr)
library(limma)
source("viz_helpers.R")

inputArgs <- commandArgs(trailingOnly=T)

dataset <- as.character(inputArgs[1])
cis.dist <- as.character(inputArgs[2])
n.samp.pcs <- as.numeric(inputArgs[3])
n.cl.pcs <- as.numeric(inputArgs[4])
cell.type.reg <- as.logical(inputArgs[5])
agg <- as.character(inputArgs[6])
chunks_tot <- as.numeric(inputArgs[7])

message(paste0("dataset, ", dataset, 
               "\ncis.dist, ", cis.dist, 
               "\nsample pcs, ", n.samp.pcs,
               "\ncell line pcs, ", n.cl.pcs,
               "\nregress props, ", cell.type.reg,
               "\naggregate, ", agg,
               "\ntotal chunks, ", chunks_tot))

stopifnot(dataset %in% c("bulk", "bulk7", "pseudobulk", "pseudobulk-cm", "pseudobulk-cf"))
stopifnot(cis.dist %in% c("10k", "25k", "50k", "100k"))
stopifnot(n.samp.pcs<=50)
stopifnot(n.cl.pcs<=10)

if (cell.type.reg) {
  cell.type.reg <- "regtypes"
} else {
  cell.type.reg <- "notypes"
}

qtl_files <- paste0("../results/eqtl_dynamic/nonlinear_dQTL/", dataset, "/", agg, "/", cis.dist, "-", n.cl.pcs, "clpcs-", n.samp.pcs, "pcs-", cell.type.reg, "-", seq(1, 1), "-t_test.tsv")
results <- map(qtl_files, read_tsv) %>% 
  bind_rows %>%
  filter(!duplicated(paste(gene, snp, sep="_")))

# Bonferroni multiple testing correction
results <- results %>%
  mutate(stdev.unscaled=stdev.unscaled.gxt2) %>%
  mutate(coefficients=beta.gxt2) %>%
  mutate(se.unadj = stdev.unscaled.gxt2*sigma) %>%
  mutate(t.unadj = beta.gxt2 / se.unadj) %>%
  mutate(p.unadj = 2*pt(-abs(t.unadj), df=df.residual[1]))

# TO REMOVE - compare to LRT method
lrt <- read_tsv("../results/eqtl_dynamic/nonlinear_dQTL/bulk/day/50k-5clpcs-0pcs-notypes-1-lrt.tsv")
comp <- results %>%
  select(gene, snp, p.unadj) %>%
  mutate(ttest=-log10(p.unadj)) %>%
  inner_join(mutate(lrt, lrtest=-log10(p)), by=c("gene", "snp"))
ggplot(comp, aes(x=ttest, y=lrtest)) + geom_point() +
  geom_abline(slope=1, intercept=0, linetype="dashed")
gene_ct <- table(results$gene)
results$bonf.p.unadj <- sapply(gene_ct[results$gene] * results$p.unadj, function(x){min(1,x)})
lrt$bonf.p.unadj <- sapply(gene_ct[results$gene] * lrt$p, function(x){min(1,x)})
top.hits.ttest  <- results %>% arrange(p.unadj) %>%
  filter(!duplicated(gene)) %>%
  mutate(qval.unadj=qvalue(bonf.p.unadj)$q)
top.hits.lrt <- lrt %>% arrange(p) %>%
  filter(!duplicated(gene)) %>%
  mutate(qval.unadj=qvalue(bonf.p.unadj)$q)

lrt.egenes <- top.hits.lrt %>% filter(qval.unadj<=0.05)
ttest.egenes <- top.hits.ttest %>% filter(qval.unadj<=0.05)

linear.egenes <- read_tsv("../results/eqtl_dynamic/linear_dQTL/bulk/day/50k-5clpcs-0pcs-notypes.tophits.tsv") %>%
  filter(qval.unadj<=0.05) %>%
  filter((gene %in% results$gene) & (snp %in% results$snp))
sum(lrt.egenes$gene %in% linear.egenes$gene)
linear.results <- read_tsv("../results/eqtl_dynamic/linear_dQTL/bulk/day/50k-5clpcs-0pcs-notypes.mtc.tsv") %>%
  filter((gene %in% results$gene) & (snp %in% results$snp))

# visualize the other LRT egenes that we are not currently detecting
lrt.unique.egenes <- setdiff(lrt.egenes$gene, c(ttest.egenes$gene, linear.egenes$gene))
egene <- lrt.unique.egenes[2]
evar <- filter(lrt.egenes, gene==!!egene)$snp
p1 <- viz_expression(egene, evar, "bulk", "day") + 
  ylab(paste0(egene))
p2 <- viz_residuals(egene, evar, "bulk", "day")  + 
  ylab(paste0(egene, " Linear Res"))
p3 <- viz_residuals_nonlinear(egene, evar, "bulk", "day")  + 
  ylab(paste0(egene, " Nonlinear Res"))
p1 + p2 + p3 + plot_layout(ncol=1)

# visualize the egenes detected by linear test that aren't detected by LRT
linear.unique.egenes <- setdiff(linear.egenes$gene, c(ttest.egenes$gene, lrt.egenes$gene))
egene <- linear.unique.egenes[2]
evar <- filter(linear.egenes, gene==!!egene)$snp
p1 <- viz_expression(egene, evar, "bulk", "day") + 
  ylab(paste0(egene))
p2 <- viz_residuals(egene, evar, "bulk", "day")  + 
  ylab(paste0(egene, " Linear Res"))
p3 <- viz_residuals_nonlinear(egene, evar, "bulk", "day")  + 
  ylab(paste0(egene, " Nonlinear Res"))
p1 + p2 + p3 + plot_layout(ncol=1)

# visualize the egenes detected by nonlinear t-test
nonlinear.t <- ttest.egenes$gene
egene <- nonlinear.t[2]
evar <- filter(ttest.egenes, gene==!!egene)$snp
p1 <- viz_expression(egene, evar, "bulk", "day") + 
  ylab(paste0(egene))
p2 <- viz_residuals(egene, evar, "bulk", "day")  + 
  ylab(paste0(egene, " Linear Res"))
p3 <- viz_residuals_nonlinear(egene, evar, "bulk", "day")  + 
  ylab(paste0(egene, " Nonlinear Res"))
p1 + p2 + p3 + plot_layout(ncol=1)

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

write_tsv(results, paste0("../results/eqtl_dynamic/linear_dQTL/", dataset, "/", agg, "/", cis.dist, "-", n.cl.pcs, "clpcs-", n.samp.pcs, "pcs-", cell.type.reg, ".mtc.tsv"))
write_tsv(top.hits, paste0("../results/eqtl_dynamic/linear_dQTL/", dataset, "/", agg, "/", cis.dist, "-", n.cl.pcs, "clpcs-", n.samp.pcs, "pcs-", cell.type.reg, ".tophits.tsv"))
