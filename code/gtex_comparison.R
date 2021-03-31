library(tidyverse)
library(vroom)

inputArgs <- commandArgs(trailingOnly=T)

dataset <- as.character(inputArgs[1])
cis.dist <- as.character(inputArgs[2])
n.samp.pcs <- as.numeric(inputArgs[3])
n.cl.pcs <- as.numeric(inputArgs[4])
cell.type.reg <- as.logical(inputArgs[5])
agg <- as.character(inputArgs[6])
qtl.type <- as.character(inputArgs[7])

# list GTEx tissues 
gtex.tissues <- list.files("../data/gtex/GTEx_Analysis_v8_eQTL") %>%
  str_extract("[^.]+") %>%
  unique

# load variant & gene dictionary 
rsid_dict <- vroom("../data/YRI.hg38.filtered.tsv", col_select=c("#CHROM", "POS", "ID")) %>%
  mutate(snp=paste(`#CHROM`, POS, sep="_"), .keep="unused") %>%
  rename(rsid=ID)
vardict <- vroom("/project2/gilad/reem/singlecellCM/dynamicqtl/gtexoverlap/variantnames.txt") %>%
  rename(rsid=`rs_id_dbSNP151_GRCh38p7`) %>%
  select(rsid, `variant_id`)
genedict <- vroom("../data/genes.tsv", col_names=c("gene_id", "gene")) 

# linear dynamic eQTLs
if (qtl.type=="linear") {
  dqtl <- vroom(paste0("../results/eqtl_dynamic/linear_dQTL/", dataset, "/", agg, "/", cis.dist, "-", n.cl.pcs, "clpcs-", n.samp.pcs, "pcs-notypes.tophits.tsv"))
} else if (qtl.type == "nonlinear") {
  dqtl <- vroom(paste0("../results/eqtl_dynamic/nonlinear_dQTL/", dataset, "/", agg, "/t_test/", cis.dist, "-", n.cl.pcs, "clpcs-", n.samp.pcs, "pcs-notypes.tophits.tsv"))
}

dqtl <- dqtl %>%
  select(gene, snp, p.unadj, bonf.p.unadj, qval.unadj) %>% 
  rename(dqtl.p=p.unadj, dqtl.p.bonf=bonf.p.unadj, dqtl.q=qval.unadj) %>%
  left_join(genedict, by="gene") %>%
  left_join(rsid_dict, by="snp") %>%
  left_join(vardict, by="rsid") %>% 
  mutate(gene_id=str_extract(gene_id, "[^.]+"))

if (qtl.type=="linear") {
  fname <- paste0("../results/eqtl_dynamic/linear_dQTL/", dataset, "/", agg, "/", cis.dist, "-", n.cl.pcs, "clpcs-", n.samp.pcs, "pcs-notypes.gtex.tsv")
} else if (qtl.type == "nonlinear") {
  fname <- paste0("../results/eqtl_dynamic/nonlinear_dQTL/", dataset, "/", agg, "/t_test/", cis.dist, "-", n.cl.pcs, "clpcs-", n.samp.pcs, "pcs-notypes.gtex.tsv")
}

if (file.exists(fname)) {
  file.remove(fname)
}

for (t in gtex.tissues) {
  gtex.overlap <- vroom(paste0("../data/gtex/GTEx_Analysis_v8_eQTL/", t, ".v8.signif_variant_gene_pairs.txt.gz")) %>%
    select(variant_id, gene_id, pval_nominal) %>%
    mutate(gene_id=str_extract(gene_id, "[^.]+")) %>%
    filter((gene_id %in% dqtl_hits$gene_id) & (variant_id %in% dqtl_hits$variant_id)) %>%
    mutate(tissue=t) %>%
    write_tsv(fname, append=T)
}

# gtex.overlap.pbcm <- vroom("../results/eqtl_dynamic/linear_dQTL/pseudobulk-cm/day/50k-5clpcs-0pcs-notypes.gtex.tsv")