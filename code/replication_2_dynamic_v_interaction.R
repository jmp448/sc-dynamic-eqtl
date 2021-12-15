library(tidyverse)

cm.ietests <- vroom("../results/eqtl_dynamic/ieQTL/bulk/CM/50k-5clpcs-0pcs-notypes.mtc.tsv") %>%
  unite(gv, gene, snp, sep="--", remove=F)
cm.ieqtls <- vroom("../results/eqtl_dynamic/ieQTL/bulk/CM/50k-5clpcs-0pcs-notypes.tophits.tsv") %>%
  unite(gv, gene, snp, sep="--", remove=F)
cm.ieqtl.hits <- filter(cm.ieqtls, qval.unadj<=0.05)

cm.tests <- vroom("../results/eqtl_dynamic/linear_dQTL/pseudobulk-cm/bin16/50k-5clpcs-0pcs-notypes.mtc.tsv") %>%
  unite(gv, gene, snp, sep="--", remove=F)
cm.dqtls <- vroom("../results/eqtl_dynamic/linear_dQTL/pseudobulk-cm/bin16/50k-5clpcs-0pcs-notypes.tophits.tsv")  %>%
  unite(gv, gene, snp, sep="--", remove=F)
cm.dqtl.hits <- filter(cm.dqtls, qval.unadj<=0.05)

length(intersect(cm.ieqtl.hits$gene, cm.dqtl.hits$gene))
cm.pi1 = 1-qvalue(filter(cm.tests, gv %in% cm.ieqtl.hits$gv)$p.unadj)$pi0
cm.pi1_rev = 1-qvalue(filter(cm.ietests, gv %in% cm.dqtl.hits$gv)$p.unadj)$pi0
cm.pi1_rev_null = 1-qvalue(sample_n(cm.ietests, 357)$p.unadj)$pi0

cf.ietests <- vroom("../results/eqtl_dynamic/ieQTL/bulk/CF/50k-5clpcs-0pcs-notypes.mtc.tsv") %>%
  unite(gv, gene, snp, sep="--", remove=F)
cf.ieqtls <- vroom("../results/eqtl_dynamic/ieQTL/bulk/CF/50k-5clpcs-0pcs-notypes.tophits.tsv") %>%
  unite(gv, gene, snp, sep="--", remove=F)
cf.ieqtl.hits <- filter(cf.ieqtls, qval.unadj<=0.05)

cf.tests <- vroom("../results/eqtl_dynamic/linear_dQTL/pseudobulk-cf/bin16/50k-5clpcs-0pcs-notypes.mtc.tsv") %>%
  unite(gv, gene, snp, sep="--", remove=F)
cf.dqtls <- vroom("../results/eqtl_dynamic/linear_dQTL/pseudobulk-cf/bin16/50k-5clpcs-0pcs-notypes.tophits.tsv") %>%
  unite(gv, gene, snp, sep="--", remove=F)
cf.dqtl.hits <- filter(cf.dqtls, qval.unadj<=0.05)

length(intersect(cf.ieqtl.hits$gene, cf.dqtl.hits$gene))
cf.pi0= 1-qvalue(filter(cf.tests, gv %in% cf.ieqtl.hits$gv)$p.unadj)$pi0
cf.pi1_rev = 1-qvalue(filter(cf.ietests, gv %in% cf.dqtl.hits$gv)$p.unadj)$pi0
cf.pi1_rev_null = 1-qvalue(sample_n(cf.ietests, 903)$p.unadj)$pi0

bulk.tests <- vroom("../results/eqtl_dynamic/linear_dQTL/bulk/day/50k-5clpcs-0pcs-notypes.mtc.tsv") %>%
  unite(gv, gene, snp, sep="--", remove=F)
cfbulk.pi1_rev = 1-qvalue(filter(bulk.tests, gv %in% cf.dqtl.hits$gv)$p.unadj)$pi0
cmbulk.pi1_rev = 1-qvalue(filter(bulk.tests, gv %in% cm.dqtl.hits$gv)$p.unadj)$pi0


