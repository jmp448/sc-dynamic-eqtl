library(tidyverse)
library(vroom)

geno <- vroom("YRI.hg38.all.vcf", skip=461) %>% 
  mutate(ref_ac=rowSums(across(starts_with("1")))) %>%
  rowwise() %>%
  mutate(mac = min(ref_ac, 38 - ref_ac)) %>%
  mutate(maf = mac / 38) %>%
  filter(maf >= 0.1) %>%
  select(!c(ref_ac, mac, maf)) %>%
  write_tsv("YRI.hg38.filtered.tsv")

geno.bed <- geno %>%
  select(`#CHROM`, POS, ID) %>%
  rename(chr=`#CHROM`, end=POS, rs_id=ID) %>%
  mutate(start=end-1, .after=1) %>%
  write_tsv("../data/genotypes/snp_locs.bed", col_names = F)

geno.freq <- vroom("../data/YRI.hg38.filtered.tsv") %>%
  mutate(ref_ac=rowSums(across(starts_with("1")))) %>%
  rowwise() %>%
  mutate(mac = min(ref_ac, 38 - ref_ac)) %>%
  ungroup() %>%
  mutate(maf = mac / 38) %>%
  mutate(snp=paste(`#CHROM`, POS, sep="_")) %>%
  select(snp, maf) %>%
  write_tsv("../data/geno_maf.tsv")


  