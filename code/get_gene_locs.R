library(tidyverse)
library(biomaRt)



# get genes and download gene locations
my_genes <- read_tsv("../data/genes.tsv", col_names=F) %>%
  `colnames<-`(c("ENSG", "HGNC"))

ensembl <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL")
ensembl <- useDataset(mart=ensembl, dataset="hsapiens_gene_ensembl")

genes <- getBM(attributes=c("ensembl_gene_id_version", "hgnc_symbol", "gene_biotype", "chromosome_name", "transcription_start_site"),
               mart=ensembl, uniqueRows = T)

# subset to protein-coding genes on somatic chromosomes
genes <- genes %>% 
  as_tibble() %>%
  filter(chromosome_name %in% as.character(seq(1,22))) %>% 
  filter(gene_biotype == "protein_coding") %>% 
  dplyr::select(-c(gene_biotype)) %>%
  rename(ENSG=ensembl_gene_id_version) %>%
  rename(HGNC=hgnc_symbol) 

my_genes <- my_genes %>%
  inner_join(genes, by=c("ENSG", "HGNC")) %>%
  filter(!duplicated(ENSG) & !duplicated(HGNC)) %>%
  write_tsv("../data/gene_locs.tsv")

my_genes_bed <- my_genes %>% 
  dplyr::select(chromosome_name, transcription_start_site, HGNC) %>%
  rename(end=transcription_start_site) %>%
  mutate(start=end-1, .after=1) %>%
  mutate(chr=paste0("chr", chromosome_name), .before=1, .keep="unused") %>%
  write_tsv("../data/gene_locs.bed", col_names=F)

