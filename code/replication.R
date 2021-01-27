library(tidyverse)
library(qvalue)

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


