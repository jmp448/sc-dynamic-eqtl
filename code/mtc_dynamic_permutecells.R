library(tidyverse)
library(stats)
library(qvalue)
library(mashr)
library(limma)
library(fitdistrplus)
library(vroom)

inputArgs <- commandArgs(trailingOnly=T)

dataset <- as.character(inputArgs[1])
cis.dist <- as.character(inputArgs[2])
n.samp.pcs <- as.numeric(inputArgs[3])
n.cl.pcs <- as.numeric(inputArgs[4])
cell.type.reg <- as.logical(inputArgs[5])
agg <- as.character(inputArgs[6])
chunks_tot <- as.numeric(inputArgs[7])
perm_samples <- as.numeric(inputArgs[8])

message(paste0("dataset, ", dataset, 
               "\ncis.dist, ", cis.dist, 
               "\nsample pcs, ", n.samp.pcs,
               "\ncell line pcs, ", n.cl.pcs,
               "\nregress props, ", cell.type.reg,
               "\naggregate, ", agg,
               "\ntotal chunks, ", chunks_tot,
               "\n# permutations, ", perm_samples))

stopifnot(dataset %in% c("bulk", "bulk7", "pseudobulk", "pseudobulk-cm", "pseudobulk-cf"))
stopifnot(cis.dist %in% c("10k", "25k", "50k", "100k"))
stopifnot(n.samp.pcs<=50)
stopifnot(n.cl.pcs<=10)

if (cell.type.reg) {
  cell.type.reg <- "regtypes"
} else {
  cell.type.reg <- "notypes"
}

load_sample <- function(s) {
  qtl_files <- paste0("../results/eqtl_dynamic/linear_dQTL/permuted_cells/", dataset, "/", agg, "/sample-", s, "/", cis.dist, "-", n.cl.pcs, "clpcs-", n.samp.pcs, "pcs-", cell.type.reg, "-", seq(1, chunks_tot), ".tsv")
  sample_tstats <- map(qtl_files, vroom) %>%
    bind_rows %>%
    filter(!duplicated(paste(gene, snp, sep="_"))) %>%
    mutate(t=beta.gxt / (stdev.unscaled.gxt * sigma)) %>%
    dplyr::select(gene, snp, t) %>%
    rename(!!paste0("sample", as.character(s)) := t)
  sample_tstats
}

get_gamma_params <- function(tstats) {
  tryCatch(
    {
      fitdist(tstats^2, "gamma")$estimate
    },
    error=function(cond){
	na.out = as.numeric(c(NA, NA))
	names(na.out) = c("shape", "rate")
	return(na.out)
    }
  )
}

merge_samples <- map(seq(1, perm_samples), load_sample) %>%
  reduce(inner_join, by=c("gene", "snp")) %>%
  write_tsv(paste0("../results/eqtl_dynamic/linear_dQTL/permuted_cells/", dataset, "/", agg, "/", cis.dist, "-", n.cl.pcs, "clpcs-", n.samp.pcs, "pcs-", cell.type.reg, ".compiled_results.tsv"))

# merge_samples <- vroom(paste0("../results/eqtl_dynamic/linear_dQTL/permuted_cells/", dataset, "/", agg, "/", cis.dist, "-", n.cl.pcs, "clpcs-", n.samp.pcs, "pcs-", cell.type.reg, ".compiled_results.tsv"))

perm_tstats <- unite(merge_samples, "gv", gene, snp, sep="--") %>%
  column_to_rownames("gv") %>%
  as.matrix
saveRDS(perm_tstats, paste0("../results/eqtl_dynamic/linear_dQTL/permuted_cells/", dataset, "/", agg, "/", cis.dist, "-", n.cl.pcs, "clpcs-", n.samp.pcs, "pcs-", cell.type.reg, ".null_tstats.tsv"))

observed_tests <- vroom(paste0("../results/eqtl_dynamic/linear_dQTL/", dataset, "/", agg, "/", cis.dist, "-", n.cl.pcs, "clpcs-", n.samp.pcs, "pcs-", cell.type.reg, ".mtc.tsv")) %>%
  dplyr::select(gene, snp, t.unadj, p.unadj) %>%
  unite("gv", gene, snp, sep="--") %>%
  filter(gv %in% rownames(perm_tstats)) %>%
  arrange(gv, rownames(perm_tstats))

emp_pvals <- empPvals(abs(observed_tests$t.unadj), abs(perm_tstats), pool=F)
emp_pvals_pooled <- empPvals(abs(observed_tests$t.unadj), abs(perm_tstats), pool=T)

# gamma_pvals <- as_tibble(t(apply(perm_tstats, 1, get_gamma_params)), rownames="gv") %>%
#   inner_join(observed_tests, by="gv") %>%
#   rename(q=t.unadj) %>%
#   dplyr::select(!p.unadj) %>%
#   column_to_rownames("gv") %>%
#   pmap(pgamma) %>%
#   unlist %>%
#   (function(x){1-x}) %>%
#   as_tibble(rownames="gv") %>%
#   write_tsv(paste0("../results/eqtl_dynamic/linear_dQTL/permuted_cells/", dataset, "/", agg, "/", cis.dist, "-", n.cl.pcs, "clpcs-", n.samp.pcs, "pcs-", cell.type.reg, ".gamma_dists.tsv"))

# results <- observed_tests %>%
#   mutate(emp.pvals=emp_pvals, gamma.pvals=gamma_pvals) %>%
#   separate(gv, c("gene", "snp"), sep="--")

results <- observed_tests %>%
  mutate(emp.pvals=emp_pvals) %>%
  mutate(emp.pvals.pooled=emp_pvals_pooled) %>%
  separate(gv, c("gene", "snp"), sep="--")

# Bonferroni multiple testing correction
gene_ct <- table(results$gene)
results$bonf.p.unadj <- sapply(gene_ct[results$gene] * results$p.unadj, function(x){min(1,x)})
results$bonf.p.emp <- sapply(gene_ct[results$gene] * results$emp.pvals, function(x){min(1,x)})
results$bonf.p.emp.pooled <- sapply(gene_ct[results$gene] * results$emp.pvals.pooled, function(x){min(1,x)})
# results$bonf.p.gamma <- sapply(gene_ct[results$gene] * results$gamma.pvals, function(x){min(1,x)})

# subset to top hits per gene and compute Storey's q value
# top.hits  <- results %>% arrange(p.unadj) %>%
#   filter(!duplicated(gene)) %>%
#   mutate(qval.unadj=qvalue(bonf.p.unadj)$q,
#          qval.emp=qvalue(bonf.p.emp)$q,
#          qval.gamma=qvalue(bonf.p.gamma)$q)
top.hits  <- results %>% arrange(p.unadj) %>%
  filter(!duplicated(gene)) %>%
  mutate(qval.unadj=qvalue(bonf.p.unadj)$q,
         qval.emp=qvalue(bonf.p.emp)$q,
         qval.emp.pooled=qvalue(bonf.p.emp.pooled)$q)

write_tsv(results, paste0("../results/eqtl_dynamic/linear_dQTL/permuted_cells/", dataset, "/", agg, "/", cis.dist, "-", n.cl.pcs, "clpcs-", n.samp.pcs, "pcs-", cell.type.reg, ".mtc.tsv"))
write_tsv(top.hits, paste0("../results/eqtl_dynamic/linear_dQTL/permuted_cells/", dataset, "/", agg, "/", cis.dist, "-", n.cl.pcs, "clpcs-", n.samp.pcs, "pcs-", cell.type.reg, ".tophits.tsv"))
