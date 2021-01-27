library(tidyverse)

hits <- read_tsv("../results/eqtl_dynamic/linear_dQTL/pseudobulk/cmbin/50k-5clpcs-0pcs.tophits.tsv")

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

test <- hits %>%
  filter(qval.unadj<=0.05) %>%
  mutate(class=pmap_chr(filter(hits, qval.unadj<=0.05), classify.dynqtl))