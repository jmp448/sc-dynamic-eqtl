library(tidyverse)
library(ashr)
library(limma)
library(stats)
library(qvalue)
library(vroom)
set.seed(1234)

inputArgs <- commandArgs(trailingOnly=T)

dataset <- as.character(inputArgs[1])
aggregate <- as.character(inputArgs[2])
cis.dist <- as.character(inputArgs[3])
npcs <- as.numeric(inputArgs[4])
cell.type <- as.character(inputArgs[5])
chunks_tot <- as.numeric(inputArgs[6])

qtl_files <- paste0("../results/eqtl_static/", dataset, "/", aggregate, "/", cell.type, "/", cis.dist, "-", npcs, "pcs-", seq(1, chunks_tot), ".tsv")
results <- map(qtl_files, vroom) %>% 
  bind_rows %>%
  filter(!duplicated(paste(gene, snp, sep="_")))

# get unadjusted p-values from a t-test
results$se.unadj <- results$sigma*results$stdev.unscaled
results$t.unadj <- results$coefficients / results$se.unadj
results$p.unadj <- 2*pt(-abs(results$t.unadj), df=results$df.residual[1])

# shrink the standard errors
results <- eBayes(results)
results$se.adj <- results$stdev.unscaled * sqrt(results$s2.post)
results <- rename(results, t.adj=t, p.adj=p.value)

# apply bonferroni correction
gene_ct <- table(results$gene)
results$bonf.p.unadj <- sapply(gene_ct[results$gene] * results$p.unadj, function(x){min(1,x)})
results$bonf.p.adj <- sapply(gene_ct[results$gene] * results$p.adj, function(x){min(1,x)})
write_tsv(results, paste0("../results/eqtl_static/", dataset, "/", aggregate, "/", cell.type, "/", cis.dist, "-", npcs, "pcs.mtc.tsv"))

# subset to top hits per gene
results <- arrange(results, p.unadj)
top.hits <- filter(results, !duplicated(results$gene))

# apply benjamini-hochberg and qvalue
top.hits$bh.bonf.unadj <- p.adjust(top.hits$bonf.p.unadj, method="fdr")
top.hits$qval.bonf.unadj <- qvalue(top.hits$bonf.p.unadj)$q

top.hits$bh.bonf.adj <- p.adjust(top.hits$bonf.p.adj, method="fdr")
top.hits$qval.bonf.adj <- qvalue(top.hits$bonf.p.adj)$q
write_tsv(top.hits, paste0("../results/eqtl_static/", dataset, "/", aggregate, "/", cell.type, "/", cis.dist, "-", npcs, "pcs.tophits.tsv"))