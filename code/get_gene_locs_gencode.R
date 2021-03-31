library(vroom)
library(tidyverse)

#! wget -O ../data/gencode.hg38.gff.gz ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.basic.annotation.gff3.gz

# load the genes we found in our dataset
my_genes <- vroom("../data/genes.tsv", col_names=c("ENSG", "HGNC")) 

# filter the annotated genome to genes we found
gencode <- vroom("../data/gencode.hg38.gff.gz", col_names=c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"), skip=7) %>%
  filter(seqname %in% paste0("chr", seq(1, 22))) %>%
  filter(feature == "gene") %>%
  dplyr::mutate(ensg=gsub(".*gene_id=([^;]+)[;].*", "\\1", attribute)) %>%
  dplyr::mutate(type=gsub(".*gene_type=([^;]+)[;].*", "\\1", attribute)) %>%
  dplyr::mutate(hgnc=gsub(".*gene_name=([^;]+)[;].*", "\\1", attribute)) %>%
  filter(type=="protein_coding") %>%
  filter(ensg %in% my_genes$ENSG)

trouble <- filter(gencode, (!ensg %in% my_genes$ENSG) & (hgnc %in% my_genes$HGNC))

# 5 HGNC symbols are duplicated here - these have overlapping loci, and will be removed from our analysis
hgnc.dup <- filter(gencode, hgnc %in% names(table(gencode$hgnc)[table(gencode$hgnc)>1]))
gencode <- gencode %>% filter(!hgnc %in% hgnc.dup$hgnc)

# save a filtered GFF file
gff.filt <- gencode %>%
  select(!c(ensg, type, hgnc)) %>%
  write_tsv("../data/gencode.hg38.filtered.gff")

# save a bed file of TSS for each gene with [0,1) indexing
bed.filt <- gencode %>%
  mutate(tss=if_else(strand=="+", start, end)) %>%
  select(seqname, tss, hgnc) %>%
  mutate(bed.start=tss-1, .after=1) %>%
  write_tsv("../data/gene_locs.bed", col_names=F)
