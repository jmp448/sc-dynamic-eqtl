library(tidyverse)
set.seed(42)

nbin <- 16
n.cl.pcs <- 5
cm_egenes <- read_tsv(paste0("../results/eqtl_dynamic/linear_dQTL/pseudobulk-cm/bin", nbin, "/50k-", n.cl.pcs, "clpcs-0pcs-notypes.tophits.tsv")) %>%
  filter(qval.unadj<=0.05) %>%
  select(gene) %>%
  write_tsv(paste0("../results/eqtl_dynamic/linear_dQTL/pseudobulk-cm/bin", nbin, "/50k-", n.cl.pcs, "clpcs-0pcs-notypes.egenes.tsv"), col_names=F)
bg_genes <- read_tsv(paste0("../results/eqtl_dynamic/linear_dQTL/pseudobulk-cm/bin", nbin, "/50k-", n.cl.pcs, "clpcs-0pcs-notypes.tophits.tsv")) %>%
  select(gene) %>%
  sample_n(nrow(cm_egenes)) %>%
  write_tsv(paste0("../results/eqtl_dynamic/linear_dQTL/pseudobulk-cm/bin", nbin, "/50k-", n.cl.pcs, "clpcs-0pcs-notypes.bg_genes.tsv"), col_names=F) 

cf_egenes <- read_tsv(paste0("../results/eqtl_dynamic/linear_dQTL/pseudobulk-cf/bin", nbin, "/50k-", n.cl.pcs, "clpcs-0pcs-notypes.tophits.tsv")) %>%
  filter(qval.unadj<=0.05) %>%
  select(gene) %>%
  write_tsv(paste0("../results/eqtl_dynamic/linear_dQTL/pseudobulk-cf/bin", nbin, "/50k-", n.cl.pcs, "clpcs-0pcs-notypes.egenes.tsv"), col_names=F)
bg_genes <- read_tsv(paste0("../results/eqtl_dynamic/linear_dQTL/pseudobulk-cf/bin", nbin, "/50k-", n.cl.pcs, "clpcs-0pcs-notypes.tophits.tsv")) %>%
  select(gene) %>%
  sample_n(nrow(cf_egenes)) %>%
  write_tsv(paste0("../results/eqtl_dynamic/linear_dQTL/pseudobulk-cf/bin", nbin, "/50k-", n.cl.pcs, "clpcs-0pcs-notypes.bg_genes.tsv"), col_names=F)

cm_enrich <- read_tsv(paste0("../results/eqtl_dynamic/linear_dQTL/pseudobulk-cm/bin16/50k-", n.cl.pcs, "clpcs-0pcs-notypes.gsea.tsv"),
                      skip=3)
cf_enrich <- read_tsv(paste0("../results/eqtl_dynamic/linear_dQTL/pseudobulk-cf/bin16/50k-", n.cl.pcs, "clpcs-0pcs-notypes.gsea.tsv"),
                      skip=3)


cm_iegenes <- read_tsv(paste0("../results/eqtl_dynamic/ieQTL/bulk/CM/50k-5clpcs-0pcs-notypes.tophits.tsv")) %>%
  filter(qval.unadj<=0.05) %>%
  select(gene) %>%
  write_tsv(paste0("../results/eqtl_dynamic/ieQTL/bulk/CM/50k-5clpcs-0pcs-notypes.egenes.tsv"), col_names=F)
bg_iegenes <- read_tsv(paste0("../results/eqtl_dynamic/ieQTL/bulk/CM/50k-5clpcs-0pcs-notypes.tophits.tsv")) %>%
  select(gene) %>%
  sample_n(nrow(cm_iegenes)) %>%
  write_tsv(paste0("../results/eqtl_dynamic/ieQTL/bulk/CM/50k-5clpcs-0pcs-notypes.bg_genes.tsv"), col_names=F)
cf_iegenes <- read_tsv(paste0("../results/eqtl_dynamic/ieQTL/bulk/CF/50k-5clpcs-0pcs-notypes.tophits.tsv")) %>%
  filter(qval.unadj<=0.05) %>%
  select(gene) %>%
  write_tsv(paste0("../results/eqtl_dynamic/ieQTL/bulk/CF/50k-5clpcs-0pcs-notypes.egenes.tsv"), col_names=F)
bg_iegenes <- read_tsv(paste0("../results/eqtl_dynamic/ieQTL/bulk/CF/50k-5clpcs-0pcs-notypes.tophits.tsv")) %>%
  select(gene) %>%
  sample_n(nrow(cf_iegenes)) %>%
  write_tsv(paste0("../results/eqtl_dynamic/ieQTL/bulk/CF/50k-5clpcs-0pcs-notypes.bg_genes.tsv"), col_names=F)
bulk_dgenes <- read_tsv(paste0("../results/eqtl_dynamic/linear_dQTL/bulk/day/50k-5clpcs-0pcs-notypes.tophits.tsv")) %>%
  filter(qval.unadj<=0.05) %>%
  select(gene) %>%
  write_tsv(paste0("../results/eqtl_dynamic/linear_dQTL/bulk/day/50k-5clpcs-0pcs-notypes.egenes.tsv"), col_names=F)
bg_dgenes <- read_tsv(paste0("../results/eqtl_dynamic/linear_dQTL/bulk/day/50k-5clpcs-0pcs-notypes.tophits.tsv")) %>%
  select(gene) %>%
  sample_n(nrow(cf_iegenes)) %>%
  write_tsv(paste0("../results/eqtl_dynamic/linear_dQTL/bulk/day/50k-5clpcs-0pcs-notypes.bg_genes.tsv"), col_names=F)
cm_ie_enrich <- read_tsv(paste0("../results/eqtl_dynamic/ieQTL/bulk/CM/50k-", n.cl.pcs, "clpcs-0pcs-notypes.gsea.tsv"),
                         skip=3)
cf_ie_enrich <- read_tsv(paste0("../results/eqtl_dynamic/ieQTL/bulk/CM/50k-", n.cl.pcs, "clpcs-0pcs-notypes.gsea.tsv"),
                      skip=3)
bulk_dqtl_enrich <- read_tsv(paste0("../results/eqtl_dynamic/linear_dQTL/bulk/day/50k-", n.cl.pcs, "clpcs-0pcs-notypes.gsea.tsv"),
                             skip=3)