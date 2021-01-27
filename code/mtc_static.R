library(tidyverse)
library(ashr)
library(limma)
library(stats)
library(qvalue)
set.seed(1234)

inputArgs <- commandArgs(trailingOnly=T)

dataset <- as.character(inputArgs[1])
aggregate <- as.character(inputArgs[2])
cis.dist <- as.character(inputArgs[3])
npcs <- as.numeric(inputArgs[4])
cell.type <- as.character(inputArgs[5])

summary.stats <- read_tsv(paste0("../results/eqtl_static/", dataset, "/", aggregate, "/", cell.type, "/", cis.dist, "-", npcs, "pcs.tsv"))

# get unadjusted p-values from a t-test
summary.stats$se.unadj <- summary.stats$sigma*summary.stats$stdev.unscaled
summary.stats$t.unadj <- summary.stats$coefficients / summary.stats$se.unadj
summary.stats$p.unadj <- 2*pt(-abs(summary.stats$t.unadj), df=summary.stats$df.residual[1])

# shrink the standard errors
summary.stats <- eBayes(summary.stats)
summary.stats$se.adj <- summary.stats$stdev.unscaled * sqrt(summary.stats$s2.post)
summary.stats <- rename(summary.stats, t.adj=t, p.adj=p.value)

# apply bonferroni correction
gene_ct <- table(summary.stats$gene)
summary.stats$bonf.p.unadj <- sapply(gene_ct[summary.stats$gene] * summary.stats$p.unadj, function(x){min(1,x)})
summary.stats$bonf.p.adj <- sapply(gene_ct[summary.stats$gene] * summary.stats$p.adj, function(x){min(1,x)})
write_tsv(summary.stats, paste0("../results/eqtl_static/", dataset, "/", aggregate, "/", cell.type, "/", cis.dist, "-", npcs, "pcs.mtc.tsv"))

# subset to top hits per gene
summary.stats <- arrange(summary.stats, p.unadj)
top.hits <- filter(summary.stats, !duplicated(summary.stats$gene))

# apply benjamini-hochberg and qvalue
top.hits$bh.bonf.unadj <- p.adjust(top.hits$bonf.p.unadj, method="fdr")
top.hits$qval.bonf.unadj <- qvalue(top.hits$bonf.p.unadj)$q

top.hits$bh.bonf.adj <- p.adjust(top.hits$bonf.p.adj, method="fdr")
top.hits$qval.bonf.adj <- qvalue(top.hits$bonf.p.adj)$q
write_tsv(top.hits, paste0("../results/eqtl_static/", dataset, "/", aggregate, "/", cell.type, "/", cis.dist, "-", npcs, "pcs.tophits.tsv"))