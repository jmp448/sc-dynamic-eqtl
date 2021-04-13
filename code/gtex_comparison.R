library(tidyverse)
library(vroom)
library(qvalue)

# load variant & gene dictionary 
rsid_dict <- vroom("../data/YRI.hg38.filtered.tsv", col_select=c("#CHROM", "POS", "ID")) %>%
  mutate(snp=paste(`#CHROM`, POS, sep="_"), .keep="unused") %>%
  rename(rsid=ID)
vardict <- vroom("/project2/gilad/reem/singlecellCM/dynamicqtl/gtexoverlap/variantnames.txt") %>%
  rename(rsid=`rs_id_dbSNP151_GRCh38p7`) %>%
  select(rsid, `variant_id`)
genedict <- vroom("../data/genes.tsv", col_names=c("gene_id", "gene")) 

# classify linear dynamic eQTLs as early/late/switch (switch included in both lists)
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

# get appropriate early/late pseudotimes (not on a 0-15 scale)
cm.medians <- read_tsv("../results/eqtl_dynamic/linear_dQTL/pseudobulk-cm/bin16/bin_medians.tsv")
t.low.cm <- cm.medians %>% 
  mutate(bin=str_extract(binind, "[^_]+$")) %>%
  filter(bin == 0) %>%
  .$t %>%
  median
t.high.cm <- cm.medians %>% 
  mutate(bin=str_extract(binind, "[^_]+$")) %>%
  filter(bin == 15) %>%
  .$t %>%
  median

cf.medians <- read_tsv("../results/eqtl_dynamic/linear_dQTL/pseudobulk-cf/bin16/bin_medians.tsv")
t.low.cf <- cf.medians %>% 
  mutate(bin=str_extract(binind, "[^_]+$")) %>%
  filter(bin == 0) %>%
  .$t %>%
  median
t.high.cf <- cf.medians %>% 
  mutate(bin=str_extract(binind, "[^_]+$")) %>%
  filter(bin == 15) %>%
  .$t %>%
  median

# 
classify.ieqtl <- function(g, intx, celltype, ..., t.low=0, t.high=1, thresh=1) {
  beta.vgt.early = (celltype*t.low + intx*0*t.low + g*0) - (celltype*t.low + intx*2*t.low + g*2)
  beta.vgt.late = (celltype*t.high + intx*0*t.high + g*0) - (celltype*t.high + intx*2*t.high + g*2)
  
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
# filter to the proper class of dynamic eQTLs (early/late/switch)
# translate rsID to variantID for comparison to GTEx, also HGNC to ENSG
convert_dqtls_linear <- function(dataset, agg, cis.dist="50k", 
                          n.cl.pcs=5, n.samp.pcs=0, 
                          class="combined", ...) {
  dqtls = vroom(paste0("../results/eqtl_dynamic/linear_dQTL/", dataset, "/", agg, "/", cis.dist, "-", n.cl.pcs, "clpcs-", n.samp.pcs, "pcs-notypes.tophits.tsv")) %>%
    filter(qval.unadj<=0.05) 
  if (class == "combined") {
    dqtls = dqtls %>%
      select(gene, snp, p.unadj, bonf.p.unadj, qval.unadj)
  } else if (class == "early") {
    if (dataset == "pseudobulk-cm") {
      dqtls = dqtls %>%
        mutate(class=pmap_chr(., classify.dynqtl, t.low=t.low.cm, t.high=t.high.cm)) %>%
        filter(class %in% c("early", "switch")) %>%
        select(gene, snp, p.unadj, bonf.p.unadj, qval.unadj)
    } else if (dataset == "pseudobulk-cf") {
      dqtls = dqtls %>%
        mutate(class=pmap_chr(., classify.dynqtl, t.low=t.low.cf, t.high=t.high.cf)) %>%
        filter(class %in% c("early", "switch")) %>%
        select(gene, snp, p.unadj, bonf.p.unadj, qval.unadj)
    } else {
      dqtls = dqtls %>%
        mutate(class=pmap_chr(., classify.dynqtl)) %>%
        filter(class %in% c("early", "switch")) %>%
        select(gene, snp, p.unadj, bonf.p.unadj, qval.unadj)
    }
  } else if (class == "late") {
    if (dataset == "pseudobulk-cm") {
      dqtls = dqtls %>%
        mutate(class=pmap_chr(., classify.dynqtl, t.low=t.low.cm, t.high=t.high.cm)) %>%
        filter(class %in% c("late", "switch")) %>%
        select(gene, snp, p.unadj, bonf.p.unadj, qval.unadj)
    } else if (dataset == "pseudobulk-cf") {
      dqtls = dqtls %>%
        mutate(class=pmap_chr(., classify.dynqtl, t.low=t.low.cf, t.high=t.high.cf)) %>%
        filter(class %in% c("late", "switch")) %>%
        select(gene, snp, p.unadj, bonf.p.unadj, qval.unadj)
    } else {
      dqtls = dqtls %>%
        mutate(class=pmap_chr(., classify.dynqtl)) %>%
        filter(class %in% c("late", "switch")) %>%
        select(gene, snp, p.unadj, bonf.p.unadj, qval.unadj)
    }
  }
  dqtls = dqtls %>%
    left_join(genedict, by="gene") %>%
    left_join(rsid_dict, by="snp") %>%
    left_join(vardict, by="rsid") %>%
    write_tsv(paste0("../data/gtex/dynamic_eqtl_results/linear_dQTL/", dataset, "/", agg, "/", class, ".tsv"))
}

# translate rsID to variantID for comparison to GTEx, also HGNC to ENSG
convert_dqtls_nonlinear <- function(dataset, agg, cis.dist="50k", 
                                    n.cl.pcs=5, n.samp.pcs=0, ...) {
  dqtls = vroom(paste0("../results/eqtl_dynamic/nonlinear_dQTL/", dataset, "/", agg, "/t_test/", cis.dist, "-", n.cl.pcs, "clpcs-", n.samp.pcs, "pcs-notypes.tophits.tsv")) %>%
    filter(qval.unadj<=0.05) %>%
    select(gene, snp, p.unadj, bonf.p.unadj, qval.unadj) %>%
    left_join(genedict, by="gene") %>%
    left_join(rsid_dict, by="snp") %>%
    left_join(vardict, by="rsid") %>%
    write_tsv(paste0("../data/gtex/dynamic_eqtl_results/nonlinear_dQTL/", dataset, "/", agg, "/combined.tsv"))
}

# now for ieQTLs
convert_ieqtls <- function(dataset, ct, cis.dist="50k", 
                                 n.cl.pcs=5, n.samp.pcs=0, 
                                 class="combined", ...) {
  ieqtls = vroom(paste0("../results/eqtl_dynamic/ieQTL/", dataset, "/", ct, "/", cis.dist, "-", n.cl.pcs, "clpcs-", n.samp.pcs, "pcs-notypes.tophits.tsv")) %>%
    filter(qval.unadj<=0.05) %>%
    rename(intx=coefficients, celltype=!!ct)
  if (class == "combined") {
    ieqtls = ieqtls %>%
      select(gene, snp, p.unadj, bonf.p.unadj, qval.unadj)
  } else if (class == "early") {
    ieqtls = ieqtls %>%
      mutate(class=pmap_chr(., classify.ieqtl)) %>%
      filter(class %in% c("early", "switch")) %>%
      select(gene, snp, p.unadj, bonf.p.unadj, qval.unadj)
  } else if (class == "late") {
    ieqtls = ieqtls %>%
      mutate(class=pmap_chr(., classify.ieqtl)) %>%
      filter(class %in% c("late", "switch")) %>%
      select(gene, snp, p.unadj, bonf.p.unadj, qval.unadj)
  }
  ieqtls = ieqtls %>%
    left_join(genedict, by="gene") %>%
    left_join(rsid_dict, by="snp") %>%
    left_join(vardict, by="rsid") %>%
    write_tsv(paste0("../data/gtex/dynamic_eqtl_results/ieQTL/", dataset, "/", ct, "/", class, ".tsv"))
}

# linear dynamic eQTLs
linear.dqtls <- tibble("dataset"=rep(c("bulk", "pseudobulk-cm", "pseudobulk-cf"), each=3), "agg"=rep(c("day", "bin16", "bin16"), each=3), "class"=rep(c("early", "late", "combined"), times=3))
pmap(linear.dqtls, convert_dqtls_linear)

# nonlinear dynamic eQTLs
nonlinear.dqtls <- tibble("dataset"=c("bulk", "pseudobulk-cm", "pseudobulk-cf"), "agg"=c("day", "bin16", "bin16"))
pmap(nonlinear.dqtls, convert_dqtls_nonlinear)

# ieQTLs
ieqtls <- tibble("dataset"=rep("bulk", 6), "ct"=rep(c("CM", "CF"), each=3), "class"=rep(c("early", "late", "combined"), times=2))
pmap(ieqtls, convert_ieqtls)

### load results and measure q-values
get_q <- function(tissue, qtl_type, dataset, agg, qtl_class, ...) {
  tryCatch(
    {vroom(paste0("../data/gtex/dynamic_eqtl_results/", qtl_type, "/", dataset, "/", agg, "/", tissue, "_", qtl_class, "_gtex_overlap.tsv")) %>%
        .$`pval_nominal` %>%
        qvalue %>%
        (function(x){1-x$pi0})},
    error=function(cond){
      return(NA)
    }
  )
}

tissues <- read_tsv("../data/gtex/dynamic_eqtl_results/tissues.txt", col_names = "tissue")
linear.dqtls <- tibble("dataset"=rep(c("bulk", "pseudobulk-cm", "pseudobulk-cf"), each=3), 
                       "agg"=rep(c("day", "bin16", "bin16"), each=3), 
                       "qtl_class"=rep(c("early", "late", "combined"), times=3),
                       "qtl_type"=rep("linear_dQTL", 9),
                       "tissue"=rep(list(tissues$tissue), 9)) %>%
  unnest(tissue) %>%
  mutate(pi1=pmap_dbl(., get_q)) %>%
  filter(linear.dqtls, qtl_class=="combined") %>%
  pivot_wider(id_cols=tissue, names_from=dataset, values_from=pi1)

nonlinear.dqtls <- tibble("dataset"=c("bulk", "pseudobulk-cm", "pseudobulk-cf"), 
                       "agg"=c("day", "bin16", "bin16"), 
                       "qtl_class"=rep("combined", 3),
                       "qtl_type"=rep("nonlinear_dQTL", 3),
                       "tissue"=rep(list(tissues$tissue), 3)) %>%
  unnest(tissue) %>%
  mutate(pi1=pmap_dbl(., get_q)) %>%
  pivot_wider(id_cols=tissue, names_from=dataset, values_from=pi1)
# note that nonlinear dynamic eQTLs have higher replication rates (in seemingly random tissues), but there are far fewer hits

ieqtls <- tibble("dataset"=rep("bulk", 6), 
                 "agg"=rep(c("CM", "CF"), each=3), 
                 "qtl_class"=rep(c("early", "late", "combined"), times=2),
                 "qtl_type"="ieQTL", 
                 "tissue"=rep(list(tissues$tissue), 6)) %>%
  unnest(tissue) %>%
  mutate(pi1=pmap_dbl(., get_q)) 

comb.ieqtl <- ieqtls %>%
  filter(qtl_class=="combined") %>%
  pivot_wider(id_cols=tissue, names_from=agg, values_from=pi1)

late.ieqtl <- ieqtls %>%
  filter(qtl_class=="late") %>%
  pivot_wider(id_cols=tissue, names_from=agg, values_from=pi1)

