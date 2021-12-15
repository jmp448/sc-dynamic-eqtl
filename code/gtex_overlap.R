library(tidyverse)

# get variants in their own column
get_var <- function(gtex_gv, ...) {
  if (is.na(gtex_gv)) {
    NA
  } else {
    s1 <- str_split(gtex_gv, "--")[[1]][2]
    chrpos <- str_split(s1, "_")[[1]][1:2]
    v <- paste(chrpos[1], chrpos[2], sep="_")
    v
  }
}

# find a set of the most transient nonlinear dynamic eQTLs
cm_nl_dqtl <- vroom("../results/eqtl_dynamic/nonlinear_dQTL/pseudobulk-cm/bin16/t_test/50k-5clpcs-0pcs-notypes.tophits.tsv") %>%
  filter(qval.unadj<0.05) %>%
  arrange(bonf.p.unadj) %>%
  unite(gv, gene, snp, sep="--", remove=F)

# beta.hat <- readRDS("../results/eqtl_static/pseudobulk/type/all_betas.rds")
# se.hat <- readRDS("../results/eqtl_static/pseudobulk/type/all_ses.rds")
# 
# beta.nldqtl <- beta.hat[rownames(beta.hat) %in% cm_nl_dqtl$gv,]
# se.nldqtl <- se.hat[rownames(se.hat) %in% cm_nl_dqtl$gv,]
# t.nldqtl <- beta.nldqtl / se.nldqtl
# 
# # get the gtex results
# linear_gtex <- read_tsv("../results/gtex_overlap/cm-linear/combined_sig_overlap_withnames.tsv") 
# linear_gtex_overlap <- linear_gtex %>%
#   mutate(gvpairs=str_split(gvpairs, ";")) %>%
#   unnest(cols=c(gvpairs)) %>%
#   drop_na
# linear_gtex_overlap$v <- sapply(linear_gtex_overlap$gvpairs, get_var)
#   
nonlinear_gtex <- read_tsv("../results/gtex_overlap/cm-nonlinear/combined_sig_overlap_withnames.tsv")
nonlinear_gtex_overlap <- nonlinear_gtex %>%
  mutate(gvpairs=str_split(gvpairs, ";")) %>%
  unnest(cols=c(gvpairs)) %>%
  drop_na
nonlinear_gtex_overlap$v <- sapply(nonlinear_gtex_overlap$gvpairs, get_var)
# 
# # how does dividing up nonlinear dynamic eQTLs into more or less transient affect the overlap?
# 
# # define transient as having lowest mean effect (t-statistic) in iPSC & CM
# endpoints <- tibble("endpoint_magnitude"=apply(abs(t.nldqtl[,c("IPSC", "CM")]), 1, mean), "gv"=rownames(t.nldqtl)) %>%
#   rowid_to_column() %>%
#   arrange(endpoint_magnitude) %>%
#   slice_head(n=34) %>%
#   mutate(v=str_extract(gv, "[^-]+$"))
# 
# # note - finding matches based on variants bc there's no eqtls that are evars for multiple genes
# transient_overlap <- filter(nonlinear_gtex_overlap, v %in% endpoints$v)
# num_overlap_transient <- length(unique(transient_overlap$v))
# nontransient_overlap <- filter(nonlinear_gtex_overlap, !v %in% endpoints$v)
# num_overlap_nontransient <- length(unique(nontransient_overlap$v))
# 
# # this isn't convincing, perhaps because we are not properly specifying transient nonlinear dynamic eQTLs
# endpoints <- tibble("endpoint_magnitude"=apply(abs(t.nldqtl[,c("CM", "CF")]), 1, mean), "gv"=rownames(t.nldqtl)) %>%
#   rowid_to_column() %>%
#   arrange(endpoint_magnitude) %>%
#   slice_head(n=34) %>%
#   mutate(v=str_extract(gv, "[^-]+$"))
# 
# transient_overlap <- filter(nonlinear_gtex_overlap, v %in% endpoints$v)
# num_overlap_transient <- length(unique(transient_overlap$v))
# nontransient_overlap <- filter(nonlinear_gtex_overlap, !v %in% endpoints$v)
# num_overlap_nontransient <- length(unique(nontransient_overlap$v))
# 
# # what if we subset GTEx tissues to only the tissues related to heart and muscle?
# heart.tissues <- c("Heart_Left_Ventricle", "Heart_Atrial_Appendage", "Artery_Aorta", "Artery_Coronary") #"Cells_Cultured_fibroblasts", "Muscle_Skeletal"
# nonlinear_gtex_overlap_heart <- filter(nonlinear_gtex_overlap, tissue %in% heart.tissues)
# linear_gtex_overlap_heart <- filter(linear_gtex_overlap, tissue %in% heart.tissues)
# 
# num_overlap_linear_heart <- length(unique(linear_gtex_overlap_heart$v))
# num_overlap_nonlinear_heart <- length(unique(nonlinear_gtex_overlap_heart$v))
# 
# transient_overlap_heart <- filter(nonlinear_gtex_overlap_heart, v %in% endpoints$v)
# num_overlap_transient_heart <- length(unique(transient_overlap_heart$v))
# nontransient_overlap_heart <- filter(nonlinear_gtex_overlap_heart, !v %in% endpoints$v)
# num_overlap_nontransient <- length(unique(nontransient_overlap$v))

# identify transient nonlinear dynamic eQTLs
classify.nl.dqtl <- function(beta.g, beta.t, beta.t2, beta.gxt, beta.gxt2, ..., t.low=0, t.middle=0.5, t.high=1, thresh=0.1) {
  beta.vgt.early = (beta.g*0 + beta.t*t.low + beta.t2*t.low^2 + beta.gxt*0*t.low + beta.gxt2 * 0*t.low^2) - (beta.g*2 + beta.t*t.low + beta.t2*t.low^2 + beta.gxt*2*t.low + beta.gxt2*2*t.low^2)
  beta.vgt.middle = (beta.g*0 + beta.t*t.middle + beta.t2*t.middle^2 + beta.gxt*0*t.middle + beta.gxt2*0*t.middle^2) - (beta.g*2 + beta.t*t.middle + beta.t2*t.middle^2 + beta.gxt*2*t.middle + beta.gxt2*2*t.middle^2)
  beta.vgt.late = (beta.g*0 + beta.t*t.high + beta.t2*t.high^2 + beta.gxt*0*t.high + beta.gxt2*0*t.high^2) - (beta.g*2 + beta.t*t.high + beta.t2*t.high^2 + beta.gxt*2*t.high + beta.gxt2*2*t.high^2)
  
  if ((abs(beta.vgt.middle)>=abs(beta.vgt.late)) & (abs(beta.vgt.middle)>=abs(beta.vgt.early))) {
    qtl.type="transient"
  } else if (sign(beta.vgt.early)==sign(beta.vgt.late)) {
    qtl.type = if_else(abs(beta.vgt.early)>=abs(beta.vgt.late), "early", "late") 
  } else if ((abs(beta.vgt.early) >= thresh) & (abs(beta.vgt.late) >= thresh)) {
      qtl.type = "switch"
  } else if ((abs(beta.vgt.early)<abs(beta.vgt.late)) & (abs(beta.vgt.early)<thresh)) {
    qtl.type = "late"
  } else if ((abs(beta.vgt.early)>abs(beta.vgt.late)) & (abs(beta.vgt.late)<thresh)) {
    qtl.type = "early"
  }
  print(paste0("early ", beta.vgt.early, ", mid ", beta.vgt.middle, ", late ", beta.vgt.late, ", type ", qtl.type))
  qtl.type
}

coefs <- vroom("../results/eqtl_dynamic/nonlinear_dQTL/pseudobulk-cm/bin16/t_test/50k-5clpcs-0pcs-notypes-coefficients.tophits.tsv")
qtl.types <- coefs %>%
  mutate(qtl.type=pmap_chr(., classify.nl.dqtl, t.low=0, t.middle=0.4, t.high=0.75, thresh=0.5)) %>%
  select(gene, snp, qtl.type) %>%
  unite(gv, gene, snp, sep="--")

cm_nl_dqtl_types <- cm_nl_dqtl %>%
  inner_join(qtl.types, by="gv")
transient_hits <- filter(cm_nl_dqtl_types, qtl.type=="transient")
early_transient_hits <- filter(cm_nl_dqtl_types, qtl.type %in% c("transient", "early"))
sum(transient_hits$snp %in% nonlinear_gtex_overlap$v)
sum(early_transient_hits$snp %in% nonlinear_gtex_overlap$v)
