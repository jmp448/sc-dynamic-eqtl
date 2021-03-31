#!/bin/bash

module load bedtools

# input - sorted bed files with (a) gene TSS, and (b) SNP locs
# output - two column bed with (1) HGNC symbol, (2) snp in chr_pos format
bedtools window -a ../data/gene_locs.sorted.bed  -b ../data/snp_locs.sorted.bed -w 50000 | \
cut -f 3,4,5,7 | \
awk 'OFS="\t"{dtss=sqrt(($4 - $1)^2); print $2, $3 "_" $4, dtss}' > ../data/all_tests.50k.d2tss.tsv