library(tidyverse)
library(qvalue)
library(vroom)

bulk.hits <- read_tsv(paste0("../data/ben/linear_dynamic_eqtls_5_pc.txt")) %>%
  filter(eFDR<=0.1) %>%
  mutate(gv=paste(rs_id, ensamble_id, sep="--"))

pseudobulk.hits <- read_tsv("/project2/gilad/reem/singlecellCM/dynamicqtl/pvaluesusingallbulktargetregions/xfinal_nF200mito30PRBDBL3_SCT_std_ensid/scdata_xfinal_dynamic_allsnpgenepairsfrombulkensid_drop700_std_celllinePCs_5covs_realdata.txt",
                            col_names=F) %>%
  rename(rs_id=X1, ensamble_id=X2, coefs=X3, p=X4) %>%
  mutate(gv=paste(rs_id, ensamble_id, sep="--"))
  
rep.rate <- pseudobulk.hits %>%
  filter(gv %in% bulk.hits$gv) %>%
  .$p %>%
  qvalue %>% 
  (function(x){1-x$pi0})

# now look just at top hit per gene
bulk.tophits <- bulk.hits %>%
  arrange(pvalue) %>%
  filter(!duplicated(ensamble_id))
rep.rate <- pseudobulk.hits %>%
  filter(gv %in% bulk.tophits$gv) %>%
  .$p %>%
  qvalue %>% 
  (function(x){1-x$pi0})

# replication of ieQTLs in bulk
n.samp.pcs <- 30
cm.ieqtl <- vroom(paste0("../results/eqtl_dynamic/ieQTL/bulk/CM/50k-5clpcs-", n.samp.pcs, "pcs-notypes.tophits.tsv")) %>%
  filter(qval.unadj<0.05)
cf.ieqtl <- vroom(paste0("../results/eqtl_dynamic/ieQTL/bulk/CF/50k-5clpcs-", n.samp.pcs, "pcs-notypes.tophits.tsv")) %>%
  filter(qval.unadj<0.05)
bulk.dynqtl <- vroom(paste0("../results/eqtl_dynamic/linear_dQTL/bulk/day/50k-5clpcs-0pcs-notypes.mtc.tsv"))

cm.pi1 <- bulk.dynqtl %>%
  filter(paste(gene, snp, sep="_") %in% paste(cm.ieqtl$gene, cm.ieqtl$snp, sep="_")) %>%
  .$p.unadj %>%
  qvalue %>%
  .$pi0 %>%
  (function(x){1-x})

cf.pi1 <- bulk.dynqtl %>%
  filter(paste(gene, snp, sep="_") %in% paste(cf.ieqtl$gene, cf.ieqtl$snp, sep="_")) %>%
  .$p.unadj %>%
  qvalue %>%
  .$pi0 %>%
  (function(x){1-x})


n.cl.pcs <- 5
cmes.ieqtl <- vroom(paste0("../results/eqtl_dynamic/ieQTL/bulk/CMES/50k-", n.cl.pcs, "clpcs-0pcs-notypes.tophits.tsv")) %>%
  filter(qval.unadj<0.05)
bulk.nl.dynqtl <- vroom(paste0("../results/eqtl_dynamic/nonlinear_dQTL/bulk/day/t_test/50k-", n.cl.pcs, "clpcs-0pcs-notypes.mtc.tsv"))
cmes.pi1 <- bulk.nl.dynqtl %>% 
  filter(paste(gene, snp, sep="_") %in% paste(cmes.ieqtl$gene, cmes.ieqtl$snp, sep="_")) %>%
  .$p.unadj %>%
  qvalue %>%
  .$pi0 %>%
  (function(x){1-x})

bulk.nl.hits <- vroom(paste0("../results/eqtl_dynamic/nonlinear_dQTL/bulk/day/t_test/50k-", n.cl.pcs, "clpcs-0pcs-notypes.tophits.tsv"))

cmes.pi2 <- bulk.nl.hits %>% 
  filter(gene %in% mes.ieqtl$gene)

n.cl.pcs <- 5
cm.dqtl <- vroom(paste0("../results/eqtl_dynamic/linear_dQTL/pseudobulk-cm/day/50k-", n.cl.pcs, "clpcs-0pcs-notypes.tophits.tsv")) %>%
  filter(qval.unadj<0.05)
cf.dqtl <- vroom(paste0("../results/eqtl_dynamic/linear_dQTL/pseudobulk-cf/day/50k-", n.cl.pcs, "clpcs-0pcs-notypes.tophits.tsv")) %>%
  filter(qval.unadj<0.05)

n.bulk.cl.pcs <- 3
bulk.dqtl <- vroom(paste0("../results/eqtl_dynamic/linear_dQTL/bulk/day/50k-", n.bulk.cl.pcs, "clpcs-0pcs-notypes.mtc.tsv"))
cm.pi1 <- bulk.dqtl %>% 
  filter(paste(gene, snp, sep="_") %in% paste(cm.dqtl$gene, cm.dqtl$snp, sep="_")) %>%
  .$p.adj %>%
  qvalue %>%
  .$pi0 %>%
  (function(x){1-x})
cf.pi1 <- bulk.dqtl %>% 
  filter(paste(gene, snp, sep="_") %in% paste(cf.dqtl$gene, cf.dqtl$snp, sep="_")) %>%
  .$p.adj %>%
  qvalue %>%
  .$pi0 %>%
  (function(x){1-x})