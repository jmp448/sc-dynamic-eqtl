library(vroom)
library(tidyverse)

classify.dynqtl <- function(beta.g, beta.gxt, beta.t, ..., t.low=0, t.high=15, thresh=1) {
  beta.vgt.early = (beta.t*t.low + beta.gxt*0*t.low + beta.g*0) - (beta.t*t.low + beta.gxt*2*t.low + beta.g*2)
  beta.vgt.late = (beta.t*t.high + beta.gxt*0*t.high + beta.g*0) - (beta.t*t.high + beta.gxt*2*t.high + beta.g*2)
  
  if (sign(beta.vgt.early)==sign(beta.vgt.late)) {
    qtl.type = if_else(abs(beta.vgt.early)>=abs(beta.vgt.late), "early", "late")
  } else {
    if ((abs(beta.vgt.early)>=abs(beta.vgt.late)) & (abs(beta.vgt.late)<thresh)) {
      qtl.type = "early"
    } else if ((abs(beta.vgt.early)<abs(beta.vgt.late)) & (abs(beta.vgt.early)<thresh)) {
      qtl.type = "late"
    } else if ((abs(beta.vgt.early)>=thresh) & (abs(beta.vgt.late)>=thresh)) {
      qtl.type = "switch"
    }
  }
}

bulk_dqtl <- vroom("../results/eqtl_dynamic/linear_dQTL/bulk7/day/50k-5clpcs-0pcs-notypes.tophits.tsv") %>%
  filter(qval.unadj<=0.05) %>%
  mutate(qtl.type=pmap_chr(., classify.dynqtl, t.low=0, t.high=15)) %>%
  select(gene, snp, qtl.type, beta.gxt, t.unadj, p.unadj, bonf.p.unadj, qval.unadj) %>%
  write_tsv("../data/dynamic-eqtl-results/bulk7_signif_linear_dqtls.tsv")

pseudobulk_dqtl <- vroom("../results/eqtl_dynamic/linear_dQTL/pseudobulk/day/50k-5clpcs-0pcs-notypes.tophits.tsv") %>%
  filter(qval.unadj<=0.05) %>%
  mutate(qtl.type=pmap_chr(., classify.dynqtl, t.low=0, t.high=15)) %>%
  select(gene, snp, qtl.type, beta.gxt, t.unadj, p.unadj, bonf.p.unadj, qval.unadj) %>%
  write_tsv("../data/dynamic-eqtl-results/pseudobulk_signif_linear_dqtls.tsv")

pseudobulk_cm_dqtl <- vroom("../results/eqtl_dynamic/linear_dQTL/pseudobulk-cm/day/50k-5clpcs-0pcs-notypes.tophits.tsv") %>%
  filter(qval.unadj<=0.05) %>%
  mutate(qtl.type=pmap_chr(., classify.dynqtl, t.low=0, t.high=15)) %>%
  select(gene, snp, qtl.type, beta.gxt, t.unadj, p.unadj, bonf.p.unadj, qval.unadj) %>%
  write_tsv("../data/dynamic-eqtl-results/pseudobulk_cm_signif_linear_dqtls.tsv")

pseudobulk_cf_dqtl <- vroom("../results/eqtl_dynamic/linear_dQTL/pseudobulk-cf/day/50k-5clpcs-0pcs-notypes.tophits.tsv") %>%
  filter(qval.unadj<=0.05) %>%
  mutate(qtl.type=pmap_chr(., classify.dynqtl, t.low=0, t.high=15)) %>%
  select(gene, snp, qtl.type, beta.gxt, t.unadj, p.unadj, bonf.p.unadj, qval.unadj) %>%
  write_tsv("../data/dynamic-eqtl-results/pseudobulk_cf_signif_linear_dqtls.tsv")