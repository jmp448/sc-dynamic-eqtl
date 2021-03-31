#!/bin/bash

wget -O ../data/gtex/eqtls.tar https://storage.googleapis.com/gtex_analysis_v8/single_tissue_qtl_data/GTEx_Analysis_v8_eQTL.tar
tar -xf ../data/gtex/eqtls.tar -C ../data/gtex

wget -O ../data/gtex/variant-dict.txt.gz https://storage.googleapis.com/gtex_analysis_v8/reference/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz