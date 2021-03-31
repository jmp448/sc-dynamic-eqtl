library(tidyverse)
library(vroom)

# load genotype info
genotypes <- vroom("../data/genotypes/YRI.hg38.filtered.tsv") %>%
  select(`#CHROM`, POS, starts_with("1")) %>%
  mutate(snp=paste(`#CHROM`, POS, sep="_"), .keep="unused", .before=1)

# filter to snps that we'll end up using
filtered_tests <- read_tsv("../data/all_tests.50k.tsv", col_names = c("gene", "snp"))

genotypes <- filter(genotypes, snp %in% filtered_tests$snp) %>%
  filter(!duplicated(snp)) %>% # no duplicate rsIDs
  write_tsv("../data/genotypes/genotypes.tsv")

# standardize genotypes
# allele_freq <- vroom("../data/genotypes/human.YRI.hg38.filtered.frq") 
# stopifnot(all.equal(allele_freq$CHROM, snp_info$chr) & all.equal(allele_freq$POS, snp_info$pos))
# allele_freq <- allele_freq %>% 
#   mutate(snp=paste(snp_info$chr, snp_info$pos, sep="_")) %>%
#   filter(snp %in% genotypes$snp) %>%
#   filter(!duplicated(snp)) %>%
#   rowwise() %>%
#   mutate(q=as.numeric(gsub(".*:(.+)\t.*", "\\1", `{ALLELE:FREQ}`))) %>%
#   mutate(p=1-q) %>%
#   select(snp, p, q)
# 
# g <- genotypes %>% select(!snp) %>% as.matrix
# std_g <- (g - 2*allele_freq$p)/sqrt(2*allele_freq$p*allele_freq$q)
# 
# standardized_genotypes <- std_g %>%
#   as_tibble() %>%
#   mutate(snp=genotypes$snp, .before=1) %>%
#   write_tsv("../data/genotypes/standardized_genotypes.tsv")
