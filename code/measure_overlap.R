library(tidyverse)
library(vroom)

all.tests <- vroom("../data/bulk/day/filtered_tests.50k.tsv")
dyn.qtls <- read_tsv("../results/eqtl_dynamic/linear_dQTL/bulk/day/50k-5clpcs-0pcs-notypes.tophits.tsv") %>%
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

combined.snps <- c(dqtl.top200$snp, matched.top200$snp)
genotypes <- filter(vroom("../data/genotypes.tsv"), snp %in% combined.snps)
write_tsv(genotypes, "temp.tsv")
# genotypes <- genotypes %>% column_to_rownames("snp") %>% t %>% 
#   as_tibble(rownames="ind")
# 
# # get a (nsnps x ncellline) matrix of each individual's genotype at each snp
# geno.dqtl <- genotypes %>% 
#   select(dqtl.top200)
# geno_comp <- function(ind1, ind2, snp, ...) {
#   geno_hits[ind1, snp] == geno_hits[ind2, snp]
# }
# tests <- bind_rows(before.top200, after.top200, before.matched, after.matched)
# geno_hits <- genotypes %>% 
#   select(ind, unique(tests$snp)) %>%
#   column_to_rownames("ind")
# inds <- rownames(geno_hits)
# tests <- tests %>%
#   slice(rep(1:n(), each=19*19)) %>%
#   mutate(ind1=rep(inds, length=nrow(.))) %>%
#   mutate(ind2=rep(inds, each=19, length=nrow(.))) %>%
#   filter(ind1!=ind2) %>%
#   mutate(same_geno=pmap_lgl(.l=., .f=geno_comp)) %>%
#   group_by(ind1, ind2, group, ct_regressed, same_geno) %>%
#   count %>%
#   filter(same_geno==T) %>%
#   mutate(pct_shared=n/200)
# 
# genotypes <- vroom("../data/genotypes.tsv") %>%
#   filter(snp %in% all.tests$snp)
