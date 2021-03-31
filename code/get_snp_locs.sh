#!/bin/bash

tail -n +466 ../data/genotypes/human.YRI.hg38.filtered.recode.vcf | \
  cut -f 1,2,3 | \
  awk 'OFS="\t"{print $1, $2-1, $2, $3}' > ../data/genotypes/snp_locs.bed