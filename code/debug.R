library(tidyverse)
library(limma)
library(vroom)
library(qvalue)
library(Seurat)

design <- read_tsv("../data/bulk/day/pcs.tsv") %>%
  mutate(ind=str_extract(sample, "[^_]+")) %>%
  mutate(day=as.numeric(str_extract(sample, "[^_]+$")))

getr2 <- function(i) {
  summary(lm(formula(paste0("PC", i, "~ind:day")), design))$adj.r.squared
}

pve <- sapply(seq(1, 50), getr2)

ggplot(tibble("pc"=as.numeric(seq(1, 50)), "r2"=pve), aes(x=pc, y=r2)) +
  geom_point() + 
  geom_segment(aes(x=pc, xend=pc, y=0, yend=r2))

# are the new hits detected after PC regression due to colliders?


get.corrs <- function(snp, pc, design) {
  summary(lm(formula(paste0(pc, "~", snp,":day")), design))$adj.r.square
}

get.comp <- function(i, pcs, geno) {
  before.hits <- read_tsv(paste0("../results/eqtl_dynamic/linear_dQTL/bulk/day/50k-5clpcs-",i-1, "pcs.tophits.tsv")) %>%
    filter(qval.unadj<=0.05) %>%
    mutate(gv=paste(gene, snp, sep="--")) %>%
    mutate(log10p=-log10(bonf.p.unadj)) %>%
    select(c(gene, snp, gv, log10p))
  
  after.hits <- read_tsv(paste0("../results/eqtl_dynamic/linear_dQTL/bulk/day/50k-5clpcs-", i, "pcs.tophits.tsv")) %>%
    filter(qval.unadj<=0.05) %>%
    mutate(gv=paste(gene, snp, sep="--")) %>%
    mutate(log10p=-log10(bonf.p.unadj)) %>%
    select(c(gene, snp, gv, log10p))
  
  diff.hits <- setdiff(after.hits$gv, before.hits$gv)
  bg.hits <- sample_n(before.hits, length(diff.hits))$gv
  
  new.hits <- filter(after.hits, gv %in% diff.hits) %>% mutate(group="new_hits")
  bg.hits <- filter(before.hits, gv %in% bg.hits) %>% mutate(group="background")
  
  all.hits <- bind_rows(new.hits, bg.hits) 
  
  u <- pcs %>% select(sample, paste0("PC", i)) %>% rename(u=starts_with("PC"))
  
  all.hits <- all.hits %>% mutate(ucorr=map_dbl(snp, get.corrs, pc=u, geno=geno)) %>% mutate(pc=i)
  
  all.hits
}

samp.pcs <- read_tsv("../data/bulk/day/pcs.tsv")
genotypes <- read_tsv("../data/genotypes.filtered.tsv") 
colnames(genotypes) <- str_replace(colnames(genotypes), "NA", "")
genotypes <- genotypes %>% column_to_rownames("snp") %>% t %>% 
  as_tibble(rownames="ind")

corrtest <- lapply(1:10, get.comp, pcs=samp.pcs, geno=genotypes) %>%
  bind_rows %>% 
  mutate(pc=factor(pc)) %>%
  write_tsv("../results/eqtl_dynamic/linear_dQTL/bulk/day/corr.analysis.tsv")

# Kclxt ~ Gxt vs u ~ Gxt
samp.pcs <- read_tsv("../data/bulk/day/pcs.tsv")
genotypes <- read_tsv("../data/genotypes.filtered.tsv") 
colnames(genotypes) <- str_replace(colnames(genotypes), "NA", "")
genotypes <- genotypes %>% column_to_rownames("snp") %>% t %>% 
  as_tibble(rownames="ind")
design <- samp.pcs %>% mutate(ind=str_extract(sample, "[^_]+")) %>%
  left_join(genotypes, by="ind")
corrtest2 <- sapply(colnames(genotypes)[2:7], get.corrs, pc=1, design=design) %>%
  as_tibble(rownames="snp") %>%
  write_tsv("../results/eqtl_dynamic/linear_dQTL/bulk/day/corr.analysis.tsv")


# get corrs and test for covariates
samp.pcs <- read_tsv("../data/bulk/day/pcs.tsv") %>% 
  mutate(ind=str_extract(sample, "[^_]+")) %>%
  mutate(day=as.numeric(str_extract(sample, "[^_]+$")))
genotypes <- read_tsv("../data/genotypes.filtered.tsv") 
colnames(genotypes) <- str_replace(colnames(genotypes), "NA", "")
genotypes <- genotypes %>% column_to_rownames("snp") %>% t %>% 
  as_tibble(rownames="ind")
pcs <- read_tsv("../data/bulk/day/clpcs.tsv") %>% 
  rename_with(function(x){str_replace(x, "PC_", "clPC")}, starts_with("PC")) %>%
  mutate(ind=as.character(ind)) %>%
  left_join(samp.pcs)
design <- pcs  %>%
  left_join(genotypes, by="ind")

clpc.str <- paste(paste0("clPC", seq(1,5)), collapse="+")

get.pve <- function(snp, pc, ..., d) {
  m <- lm(formula(paste0(pc, "~(", clpc.str, "+", snp, ")*day")), d)
  as_tibble(anova(m), rownames="coef") %>%
    filter(coef != "Residuals") %>%
    mutate(pve=100*`Sum Sq`/sum(`Sum Sq`)) %>%
    filter(coef==paste0(snp, ":day")) %>%
    .$pve
}

before.hits <- read_tsv(paste0("../results/eqtl_dynamic/linear_dQTL/bulk/day/50k-0clpcs-0pcs.tophits.tsv")) %>%
  filter(qval.unadj<=0.05) %>%
  mutate(gv=paste(gene, snp, sep="--")) %>%
  mutate(log10p=-log10(bonf.p.unadj)) %>%
  select(c(gene, snp, gv, log10p))

after.hits <- read_tsv(paste0("../results/eqtl_dynamic/linear_dQTL/bulk/day/50k-0clpcs-20pcs.tophits.tsv")) %>%
  filter(qval.unadj<=0.05) %>%
  mutate(gv=paste(gene, snp, sep="--")) %>%
  mutate(log10p=-log10(bonf.p.unadj)) %>%
  select(c(gene, snp, gv, log10p))

diff.hits <- setdiff(after.hits$gv, before.hits$gv) 
bg.hits <- sample_n(before.hits, min(length(diff.hits), nrow(before.hits)))$gv

numpcs <- 20
tests <- bind_rows(tibble("gv"=diff.hits, "group"="new_hits"),
                   tibble("gv"=bg.hits, "group"="background")) %>%
  mutate(snp=str_replace(gv, ".*--", "")) %>%
  group_by(group) %>%
  slice(rep(1:n(), each=numpcs)) %>%
  ungroup() %>%
  mutate(pc=rep(paste0("PC", seq(1:numpcs)), length=nrow(.))) %>%
  mutate(pc=factor(pc, levels=paste0("PC", seq(1:numpcs))))
design_filt <- select(design, ind, sample, day,
                      starts_with("PC"), starts_with("clPC"),
                      unique(tests$snp))
tests <- mutate(tests, pve=pmap_dbl(.l=tests, .f=get.pve, d=design_filt)) %>%
  mutate(pc=factor(pc, levels=paste0("PC", seq(1:numpcs)))) %>%
  write_tsv("../results/eqtl_dynamic/linear_dQTL/bulk/day/corr.analysis.0cl-20pcs_regressed.tsv")

ggplot(tests, aes(x=pc,y=pve,fill=group)) +
  geom_boxplot()

tests2 <- slice_head(tests, n=20)
system.time(mutate(tests2, pve=pmap_dbl(.l=tests2, .f=get.pve, d=design_filt)))

### look at how many cell lines stick together
genotypes <- vroom("../data/genotypes.tsv")
genotypes <- genotypes %>% column_to_rownames("snp") %>% t %>% 
  as_tibble(rownames="ind")

before.hits <- read_tsv("../results/eqtl_dynamic/linear_dQTL/bulk/day/50k-5clpcs-0pcs-notypes.tophits.tsv") %>%
  filter(qval.unadj<=0.05) %>%
  mutate(gv=paste(gene, snp, sep="--"))
after.hits <- read_tsv("../results/eqtl_dynamic/linear_dQTL/bulk/day/50k-5clpcs-0pcs-regtypes.tophits.tsv") %>%
  filter(lfsr<=0.05) %>%
  mutate(gv=paste(gene, snp, sep="--"))

all.tests <- read_tsv("../data/gv_pairs.filtered.50k.tsv") %>%
  mutate(gv=paste(gene, snp, sep="--"))
get.matched <- function(tests) {
  elim.genes <- str_extract(tests, "[^-]+")
  match.dists <- all.tests %>% 
    filter(gv %in% tests) %>%
    .$dist2tss
  candidates <- all.tests %>%
    filter(!gene %in% elim.genes)
  matched.tests <- bind_rows(map_dfr(match.dists, function(d){slice(candidates, which.min(abs(dist2tss-d)))}))$gv
  matched.tests
}


before.top200 <- before.hits %>% 
  arrange(lfsr) %>% 
  slice_head(n=200) %>% 
  select(gene, snp, gv) %>%
  mutate(group="dynqtl") %>%
  mutate("ct_regressed"=F)
after.top200 <- after.hits %>% 
  filter(!gv %in% before.hits$gv) %>%
  arrange(lfsr) %>% 
  slice_head(n=200) %>% 
  select(gene, snp, gv) %>% 
  mutate(group="dynqtl") %>%
  mutate("ct_regressed"=T)
before.matched <- get.matched(before.top200$gv) %>%
  as_tibble() %>%
  rename(gv=value) %>%
  mutate(gene=str_extract(gv, "[^-]+")) %>%
  mutate(snp=str_extract(gv, "[^-]+$")) %>%
  mutate(group="background") %>%
  mutate("ct_regressed"=F)
after.matched <- get.matched(after.top200$gv) %>%
  as_tibble() %>%
  rename(gv=value) %>%
  mutate(gene=str_extract(gv, "[^-]+")) %>%
  mutate(snp=str_extract(gv, "[^-]+$")) %>%
  mutate(group="background") %>%
  mutate("ct_regressed"=T)

geno_comp <- function(ind1, ind2, snp, ...) {
  geno_hits[ind1, snp] == geno_hits[ind2, snp]
}
tests <- bind_rows(before.top200, after.top200, before.matched, after.matched)
geno_hits <- genotypes %>% 
  select(ind, unique(tests$snp)) %>%
  column_to_rownames("ind")
inds <- rownames(geno_hits)
tests <- tests %>%
  slice(rep(1:n(), each=19*19)) %>%
  mutate(ind1=rep(inds, length=nrow(.))) %>%
  mutate(ind2=rep(inds, each=19, length=nrow(.))) %>%
  filter(ind1!=ind2) %>%
  mutate(same_geno=pmap_lgl(.l=., .f=geno_comp)) %>%
  group_by(ind1, ind2, group, ct_regressed, same_geno) %>%
  count %>%
  filter(same_geno==T) %>%
  mutate(pct_shared=n/200)
  
# remove duplicate rows
tests_nodub <- tests %>%
  mutate(ind1=as.integer(ind1)) %>%
  mutate(ind2=as.integer(ind2)) %>%
  mutate(pair=paste(min(ind1, ind2), "--", max(ind1, ind2))) %>%
  group_by(group, ct_regressed, same_geno) %>%
  filter(!duplicated(pair))

ggplot(tests_nodub, aes(x=ct_regressed, y=pct_shared, fill=group)) +
  geom_violin()

# paired t test - first, are differences normally distributed
reg_replace <- function(ct_regressed,...) {
  if (ct_regressed) {
    "with_regression"
  } else {
    "without_regression"
  }
}
t_testin <- tests_nodub %>%
  ungroup %>%
  mutate(reg=pmap_chr(.l=., .f=reg_replace)) %>%
  select(!c(n, ct_regressed)) %>%
  spread(key=reg, value=pct_shared) %>%
  mutate(diff=with_regression-without_regression)
ggplot(t_testin, aes(x=diff)) +
  geom_density() # not too normal, but should be ok
t.test(t_testin$without_regression, t_testin$with_regression, paired=T, alternative="greater")

# look at heat map
hm.order <- tests %>% filter(group=="dynqtl", ct_regressed==T, same_geno==T) %>%
  ungroup %>%
  select(ind1, ind2, pct_shared) %>%
  spread(key=ind2, value=pct_shared) %>%
  column_to_rownames("ind1") %>%
  dist %>%
  hclust %>%
  .$order
testit <- tests %>%
  mutate(ind1=factor(as.character(ind1), levels=ind1[hm.order])) %>%
  mutate(ind2=factor(as.character(ind2), levels=levels(ind1)))
ggplot(testit, aes(x=ind1, y=ind2, fill=pct_shared)) +
  geom_tile() + 
  facet_grid(rows=vars(group), cols=vars(ct_regressed))


my_old_genotypes <- read_tsv("../data/gv_pairs.filtered.25k.tsv")

ben_results <- read_tsv("../data/ben/linear_dynamic_eqtls_5_pc.txt")

my_results <- vroom("../data/filtered_tests.50k.tsv")

snp_translator <- vroom("../data/genotypes/snp_locs.bed", col_names=c("chr", "start", "stop", "rs_id")) %>%
  dplyr::select(!start) %>%
  mutate(snp=paste(chr, stop, sep="_"), .keep="unused")
gene_translator <- read_tsv("../data/genes.tsv", col_names = c("ensamble_id", "gene")) %>%
  dplyr::mutate(ensamble_id=str_extract(`ensamble_id`, "[^.]+"))

my_results <- left_join(my_results, snp_translator, by="snp") %>%
  left_join(gene_translator, by="gene")

shared.genes <- intersect(my_results$ensamble_id, ben_results$ensamble_id)
me_shared <- filter(my_results, ensamble_id %in% shared.genes)
ben_shared <- filter(ben_results, ensamble_id %in% shared.genes)

mutual <- inner_join(ben_shared, me_shared, by=c("rs_id", "ensamble_id"))
ben_unique <- filter(ben_shared, !rs_id %in% me_shared$rs_id)
me_unique <- filter(me_shared, !rs_id %in% ben_shared$rs_id)

sharing_dist <- ben_shared %>% mutate(snp_rep=rs_id %in% me_shared$rs_id) %>%
  group_by(ensamble_id) %>%
  summarize(shared=sum(snp_rep)/table(ensamble_id))

sharing_dist <- me_shared %>% mutate(snp_rep=rs_id %in% ben_shared$rs_id) %>%
  group_by(ensamble_id) %>%
  summarize(shared=sum(snp_rep)/table(ensamble_id))

ggplot(sharing_dist, aes(x=shared)) + geom_density()

bulk_permutations <- read_tsv("../results/eqtl_dynamic/linear_dQTL/bulk/day/50k-5clpcs-0pcs.permuted.mtc.tsv") %>%
  dplyr::select(p.unadj) %>%
  arrange(p.unadj) %>%
  rowid_to_column() %>%
  mutate(unif=rowid/nrow(.)) %>%
  mutate(log10p=-log10(p.unadj)) %>%
  mutate(log10unif=-log10(unif))

ggplot(permutations, aes(x=log10unif, y=log10p)) + 
  geom_point() +
  geom_abline(aes(slope=1, intercept=0, colour="red", linetype="dashed"))

i=6
chunk_size <- read_tsv(paste0("../results/eqtl_dynamic/ieQTL/pseudobulk/CM/50k-5clpcs-0pcs-notypes-", i, ".tsv")) %>%
  nrow

# how much variation in pseudotime is explained by individual effects? 03/09/2021
# single cell level
sc_cf <- readRDS("../data/seurat.cf.rds")
meta_cf <- as_tibble(sc_cf@meta.data, rownames="cell")
model_cf <- lm(pseudotime ~ individual, meta_cf)
summary(model_cf)$adj.r.squared

sc_cm <- readRDS("../data/seurat.cm.rds")
meta_cm <- as_tibble(sc_cm@meta.data, rownames="cell")
model_cm <- lm(pseudotime ~ individual, meta_cm)

# after aggregation
sc_cf <- read_tsv(paste0("../results/eqtl_dynamic/linear_dQTL/pseudobulk-cf/bin15/bin_medians.tsv")) %>%
  rename(sample=binind) %>%
  mutate(ind=str_sub(sample, 1, 5))
model_cf <- lm(t ~ ind, sc_cf)
summary(model_cf)$adj.r.squared

sc_cm <- read_tsv(paste0("../results/eqtl_dynamic/linear_dQTL/pseudobulk-cm/bin15/bin_medians.tsv")) %>%
  rename(sample=binind) %>%
  mutate(ind=str_sub(sample, 1, 5))
model_cm <- lm(t ~ ind, sc_cm)
summary(model_cm)$adj.r.squared


# cell line PC trends
analyses <- tibble("dataset"=rep(c("bulk", "pseudobulk-cm", "pseudobulk-cf"), each=5),
                   "agg"=c(rep("day", 5), rep("bin16", 10)),
                   "n.cl.pcs"=rep(seq(1,5), times=3))

count.hits <- function(dataset, agg, n.cl.pcs, ...) {
  nhits <- vroom(paste0("../results/eqtl_dynamic/linear_dQTL/", dataset, "/", agg, "/50k-", n.cl.pcs, "clpcs-0pcs-notypes.tophits.tsv")) %>%
    filter(qval.unadj<=0.05) %>%
    nrow(.)
  nhits
}

analyses <- analyses %>% mutate(nhits=pmap_dbl(.l=., .f=count.hits))

ggplot(analyses, aes(x=n.cl.pcs, y=nhits)) +
  geom_bar(stat="identity") +
  facet_grid(rows=vars(dataset))

egenes.5 <- vroom(paste0("../results/eqtl_dynamic/linear_dQTL/pseudobulk-cf/bin16/50k-5clpcs-0pcs-notypes.tophits.tsv")) %>%
  filter(qval.unadj<=0.05)
egenes.1 <- vroom(paste0("../results/eqtl_dynamic/linear_dQTL/pseudobulk-cf/bin16/50k-1clpcs-0pcs-notypes.tophits.tsv")) %>%
  filter(qval.unadj<=0.05)

sum(egenes.5$gene %in% egenes.1$gene)

results.5 <- vroom(paste0("../results/eqtl_dynamic/linear_dQTL/pseudobulk-cf/bin16/50k-5clpcs-0pcs-notypes.mtc.tsv")) %>%
  filter((gene %in% egenes.5$gene) & (snp %in% egenes.5$snp))
results.1 <- vroom(paste0("../results/eqtl_dynamic/linear_dQTL/pseudobulk-cf/bin16/50k-1clpcs-0pcs-notypes.mtc.tsv")) %>%
  filter((gene %in% egenes.5$gene) & (snp %in% egenes.5$snp))

p.comp <- inner_join(select(mutate(results.5, cl5=-log10(p.unadj)), c(gene, snp, cl5)),
                     select(mutate(results.1, cl1=-log10(p.unadj)), c(gene, snp, cl1)),
                     by=c("gene", "snp")) %>%
  mutate(significant.1=gene %in% egenes.1$gene)

ggplot(p.comp, aes(x=cl1, y=cl5, color=significant.1)) + 
  geom_point() + 
  geom_abline(intercept=0, slope=1, linetype="dashed") +
  xlab("-log10(P), 1 CLPC") + 
  ylab("-log10(P), 5 CLPCs") +
  labs(color="Significant, 1 CLPC") +
  theme_classic()

### look at how many cell lines stick together - LATE VERSION
dataset <- "pseudobulk-cm"
agg <- "bin16"
all.tests <- vroom(paste0("../data/", dataset, "/", agg, "/filtered_tests.50k.tsv"))
dyn.qtls <- read_tsv(paste0("../results/eqtl_dynamic/linear_dQTL/", dataset, "/", agg, "/50k-5clpcs-0pcs-notypes.tophits.tsv")) %>%
  filter(qval.unadj<=0.05) %>%
  mutate(gv=paste(gene, snp, sep="--"))

mafs <- vroom("../data/geno_maf.tsv")
dqtl.top200 <- dyn.qtls %>% 
  arrange(qval.unadj) %>% 
  slice_head(n=200) %>% 
  select(gene, snp, gv) %>%
  mutate(group="dynqtl") %>%
  mutate("ct_regressed"=F) %>%
  left_join(mafs, by="snp") %>%
  mutate(maf.bin=floor(maf/0.05))
bin.counts <- bin.counts <- count(dqtl.top200, maf.bin)
dqtl.top200 <- select(dqtl.top200, snp) 
# match MAF for some background variants
matched.top200 <- all.tests %>%
  left_join(mafs, by="snp")%>%
  mutate(maf.bin=floor(maf/0.05)) %>%
  group_by(maf.bin) %>%
  nest() %>%
  ungroup() %>%
  left_join(bin.counts, by="maf.bin") %>%
  replace_na(list(n=0)) %>%
  mutate(samp=map2(data, n, sample_n)) %>%
  select(!data) %>%
  unnest(samp) %>%
  select(snp)

# get a (nsnps x ncellline) matrix of each individual's genotype at each snp
genotypes <- vroom("../data/genotypes.tsv")
dqtl.cor <- genotypes %>% 
  filter(snp %in% dqtl.top200$snp) %>%
  column_to_rownames("snp") %>%
  cor
matched.cor <- genotypes %>%
  filter(snp %in% matched.top200$snp) %>%
  column_to_rownames("snp") %>%
  cor 

corr.cols <- colorRampPalette(brewer.pal(11, "RdBu"))(25)

png(paste0('../figs/supp/dqtl_corr_', dataset, '_', agg, '_hits.png'), width=800, height=600)
heatmap.2(dqtl.cor, Rowv=T, Colv=T, dendrogram='none', trace='none', col=rev(corr.cols), breaks=seq(-1, 1, length.out=26))
dev.off()

png(paste0('../figs/supp/dqtl_corr_', dataset, '_', agg, '_matched.png'), width=800, height=600)
heatmap.2(matched.cor, Rowv=T, Colv=T, dendrogram='none', trace='none', col=rev(corr.cols), breaks=seq(-1, 1, length.out=26))
dev.off()

#old stuff
geno_comp <- function(ind1, ind2, snp, ...) {
  geno_hits[ind1, snp] == geno_hits[ind2, snp]
}
tests <- bind_rows(before.top200, after.top200, before.matched, after.matched)
geno_hits <- genotypes %>% 
  select(ind, unique(tests$snp)) %>%
  column_to_rownames("ind")
inds <- rownames(geno_hits)
tests <- tests %>%
  slice(rep(1:n(), each=19*19)) %>%
  mutate(ind1=rep(inds, length=nrow(.))) %>%
  mutate(ind2=rep(inds, each=19, length=nrow(.))) %>%
  filter(ind1!=ind2) %>%
  mutate(same_geno=pmap_lgl(.l=., .f=geno_comp)) %>%
  group_by(ind1, ind2, group, ct_regressed, same_geno) %>%
  count %>%
  filter(same_geno==T) %>%
  mutate(pct_shared=n/200)

# remove duplicate rows
tests_nodub <- tests %>%
  mutate(ind1=as.integer(ind1)) %>%
  mutate(ind2=as.integer(ind2)) %>%
  mutate(pair=paste(min(ind1, ind2), "--", max(ind1, ind2))) %>%
  group_by(group, ct_regressed, same_geno) %>%
  filter(!duplicated(pair))

ggplot(tests_nodub, aes(x=ct_regressed, y=pct_shared, fill=group)) +
  geom_violin()


### are clpc interaction collider effects a problem? 
# for pcs 1, 2, 3, 4, 5, are new hits detected after cell line PC regression correlated w the intx term?
clpc <- 5
# get clpc * time
colliders <- vroom("../results/eqtl_dynamic/linear_dQTL/pseudobulk-cf/bin16/bin_medians.tsv") %>%
  mutate(ind=str_sub(binind, 1, 5)) %>%
  inner_join(mutate(vroom("../data/pseudobulk-cf/bin16/clpcs.tsv"), ind=as.character(ind)), by="ind")
clxt <- select(colliders, c(ind, t, paste0("PC_", clpc))) %>%
  mutate(clxt=t*!!as.name(paste0("PC_", clpc))) 
# get snps before and after pc regression
mafs <- vroom("../data/geno_maf.tsv")
hits.pre <- vroom(paste0("../results/eqtl_dynamic/linear_dQTL/pseudobulk-cf/bin16/50k-", clpc-1, "clpcs-0pcs-notypes.tophits.tsv")) %>%
  filter(qval.unadj<=0.05) %>%
  left_join(mafs, by="snp") %>%
  mutate(maf.bin=floor(maf/0.05))
hits.post <- vroom(paste0("../results/eqtl_dynamic/linear_dQTL/pseudobulk-cf/bin16/50k-", clpc, "clpcs-0pcs-notypes.tophits.tsv")) %>%
  filter(qval.unadj<=0.05) %>%
  filter(!gene %in% hits.pre$gene) %>%
  sample_n(min(nrow(.), min(100, nrow(hits.pre)))) %>%
  left_join(mafs, by="snp") %>%
  mutate(maf.bin=floor(maf/0.05))
bin.counts <- count(hits.post, maf.bin)
hits.post <- hits.post %>%
  select(snp) %>%
  mutate(grp="post")
# match pre-regression hits to an MAF bin
hits.pre <- hits.pre %>%
  group_by(maf.bin) %>%
  nest() %>%
  ungroup() %>%
  left_join(bin.counts, by="maf.bin") %>%
  replace_na(list(n=0)) %>%
  mutate(samp=map2(data, n, sample_n)) %>%
  select(!data) %>%
  unnest(samp) %>%
  select(snp) %>%
  mutate(grp="pre")
snps <- bind_rows(hits.pre, hits.post)
# load genotypes for the relevant snps
collider.check <- vroom("../data/genotypes.tsv") %>%
  filter(snp %in% snps$snp) %>%
  column_to_rownames("snp") %>% t %>% as_tibble(rownames="ind") %>%
  right_join(clxt, by="ind") %>% select(starts_with("chr"))
snps$collider.corr <- abs(apply(collider.check, 2, function(x){cor(clxt$clxt, x)}))[snps$snp]
snps$grp=factor(snps$grp, levels=c("pre", "post"))

ggplot(snps, aes(x=grp, y=collider.corr, fill=grp)) + 
  geom_boxplot() +
  theme_classic()
t.test(filter(snps, grp=="pre")$collider.corr, filter(snps, grp=="post")$collider.corr, alternative="less")
  

### look at this same issue in bulk
# for pcs 2, 3, 4, 5, are new hits detected after cell line PC regression correlated w the intx term?
clpc <- 5
# get clpc * time
colliders <- vroom("../data/bulk/day/pcs.tsv") %>%
  select(sample) %>%
  mutate(ind=str_sub(sample, 1, 5)) %>%
  mutate(t=as.numeric(str_extract(sample, "[^_]+$"))) %>%
  left_join(mutate(vroom("../data/bulk/day/clpcs.tsv"), ind=as.character(ind)), by="ind")
clxt <- select(colliders, c(ind, t, paste0("PC_", clpc))) %>%
  mutate(clxt=t*!!as.name(paste0("PC_", clpc))) 
# get snps before and after pc regression and assign to MAF bins
hits.pre <- vroom(paste0("../results/eqtl_dynamic/linear_dQTL/bulk/day/50k-", clpc-1, "clpcs-0pcs-notypes.tophits.tsv")) %>%
  filter(qval.unadj<=0.05) %>%
  left_join(mafs, by="snp") %>%
  mutate(maf.bin=floor(maf/0.05))
hits.post <- vroom(paste0("../results/eqtl_dynamic/linear_dQTL/bulk/day/50k-", clpc, "clpcs-0pcs-notypes.tophits.tsv")) %>%
  filter(qval.unadj<=0.05) %>%
  left_join(mafs, by="snp") %>%
  mutate(maf.bin=floor(maf/0.05)) %>%
  filter(!gene %in% hits.pre$gene) %>%
  sample_n(min(nrow(.), min(100, nrow(hits.pre)))) 
bin.counts <- count(hits.post, maf.bin)
hits.post <- hits.post %>%
  select(snp) %>%
  mutate(grp="post")
# match pre-regression hits to an MAF bin
hits.pre <- hits.pre %>%
  group_by(maf.bin) %>%
  nest() %>%
  ungroup() %>%
  left_join(bin.counts, by="maf.bin") %>%
  replace_na(list(n=0)) %>%
  mutate(samp=map2(data, n, sample_n)) %>%
  select(!data) %>%
  unnest(samp) %>%
  select(snp) %>%
  mutate(grp="pre")
snps <- bind_rows(hits.pre, hits.post)
# load genotypes for the relevant snps
collider.check <- vroom("../data/genotypes.tsv") %>%
  filter(snp %in% snps$snp) %>%
  column_to_rownames("snp") %>% t %>% as_tibble(rownames="ind") %>%
  right_join(clxt, by="ind") %>% select(starts_with("chr"))
snps$collider.corr <- abs(apply(collider.check, 2, function(x){cor(clxt$clxt, x)}))[snps$snp]
snps$grp=factor(snps$grp, levels=c("pre", "post"))

ggplot(snps, aes(x=grp, y=collider.corr, fill=grp)) + 
  geom_boxplot() +
  theme_classic()

t.test(filter(snps, grp=="pre")$collider.corr, filter(snps, grp=="post")$collider.corr, alternative="less")

# same issue, nonlinear dynamic eQTLs
clpc <- 1
# get clpc * time
colliders <- vroom("../results/eqtl_dynamic/linear_dQTL/pseudobulk-cf/bin16/bin_medians.tsv") %>%
  mutate(ind=str_sub(binind, 1, 5)) %>%
  inner_join(mutate(vroom("../data/pseudobulk-cf/bin16/clpcs.tsv"), ind=as.character(ind)), by="ind")
clxt <- select(colliders, c(ind, t, paste0("PC_", clpc))) %>%
  mutate(clxt=t*!!as.name(paste0("PC_", clpc))) 
# get snps before and after pc regression
mafs <- vroom("../data/geno_maf.tsv")
hits.pre <- vroom(paste0("../results/eqtl_dynamic/nonlinear_dQTL/pseudobulk-cf/bin16/t_test/50k-", clpc-1, "clpcs-0pcs-notypes.tophits.tsv")) %>%
  filter(qval.unadj<=0.05) %>%
  left_join(mafs, by="snp") %>%
  mutate(maf.bin=floor(maf/0.05))
hits.post <- vroom(paste0("../results/eqtl_dynamic/nonlinear_dQTL/pseudobulk-cf/bin16/t_test/50k-", clpc, "clpcs-0pcs-notypes.tophits.tsv")) %>%
  filter(qval.unadj<=0.05) %>%
  filter(!gene %in% hits.pre$gene) %>%
  sample_n(min(nrow(.), min(100, nrow(hits.pre)))) %>%
  left_join(mafs, by="snp") %>%
  mutate(maf.bin=floor(maf/0.05))
bin.counts <- count(hits.post, maf.bin)
hits.post <- hits.post %>%
  select(snp) %>%
  mutate(grp="post")
# match pre-regression hits to an MAF bin
hits.pre <- hits.pre %>%
  group_by(maf.bin) %>%
  nest() %>%
  ungroup() %>%
  left_join(bin.counts, by="maf.bin") %>%
  replace_na(list(n=0)) %>%
  mutate(samp=map2(data, n, sample_n)) %>%
  select(!data) %>%
  unnest(samp) %>%
  select(snp) %>%
  mutate(grp="pre")
snps <- bind_rows(hits.pre, hits.post)
# load genotypes for the relevant snps
collider.check <- vroom("../data/genotypes.tsv") %>%
  filter(snp %in% snps$snp) %>%
  column_to_rownames("snp") %>% t %>% as_tibble(rownames="ind") %>%
  right_join(clxt, by="ind") %>% select(starts_with("chr"))
snps$collider.corr <- abs(apply(collider.check, 2, function(x){cor(clxt$clxt, x)}))[snps$snp]
snps$grp=factor(snps$grp, levels=c("pre", "post"))

ggplot(snps, aes(x=grp, y=collider.corr, fill=grp)) + 
  geom_boxplot() +
  theme_classic()
t.test(filter(snps, grp=="pre")$collider.corr, filter(snps, grp=="post")$collider.corr, alternative="less")


# interaction eqtls vs bulk dynamic eqtls
bulk.linear.dqtl <- vroom("../results/eqtl_dynamic/linear_dQTL/bulk/day/50k-3clpcs-0pcs-notypes.tophits.tsv")
bulk.linear.dqtl.egenes <- filter(bulk.linear.dqtl, qval.unadj<=0.05)$gene
bulk.nonlinear.dqtl <- vroom("../results/eqtl_dynamic/nonlinear_dQTL/bulk/day/t_test/50k-5clpcs-0pcs-notypes.tophits.tsv")
bulk.nonlinear.dqtl.egenes <- filter(bulk.nonlinear.dqtl, qval.unadj<=0.05)$gene
bulk.cm.ieqtl <- vroom("../results/eqtl_dynamic/ieQTL/bulk/CM/50k-5clpcs-0pcs-notypes.tophits.tsv")
bulk.cm.ieqtl.egenes <- filter(bulk.cm.ieqtl, qval.unadj<=0.05)$gene
bulk.cf.ieqtl <- vroom("../results/eqtl_dynamic/ieQTL/bulk/CF/50k-5clpcs-0pcs-notypes.tophits.tsv")
bulk.cf.ieqtl.egenes <- filter(bulk.cf.ieqtl, qval.unadj<=0.05)$gene


sc <- readRDS("../data/seurat.normalized.rds")

# check gtex overlap 
reem.overlaps <- read_tsv("/project2/gilad/reem/singlecellCM/dynamicqtl/gtexoverlap/output_eqtloverlap/")
