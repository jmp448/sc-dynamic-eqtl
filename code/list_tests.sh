#!/bin/bash

module load bedtools

bedtools sort -i ../data/genotypes/snp_locs.bed > ../data/snp_locs.sorted.bed
bedtools sort -i ../data/gene_locs.bed > ../data/gene_locs.sorted.bed

# input - sorted bed files with (a) gene TSS, and (b) SNP locs
# output - two column bed with (1) HGNC symbol, (2) snp in chr_pos format
bedtools window -a ../data/gene_locs.sorted.bed  -b ../data/snp_locs.sorted.bed -w 50000 | \
cut -f 4,5,7 | \
awk 'OFS="\t"{print $1, $2 "_" $3}' > ../data/all_tests.50k.tsv