#!/bin/bash
source activate /project2/gilad/jpopp/sc-dynamic-eqtl/code/gsea/py_gsea
# define the gene sets
hallmark="genesets/h.all.v7.0.symbols.gmt.txt"

# agg="bin16"
# # for dataset in "pseudobulk-cm" "pseudobulk-cf" "pseudobulk"
# for dataset in "pseudobulk-cm" "pseudobulk-cf"
# do
#   test="../../results/eqtl_dynamic/linear_dQTL/${dataset}/${agg}/50k-3clpcs-0pcs-notypes.egenes.tsv"
#   background="../../results/eqtl_dynamic/linear_dQTL/${dataset}/${agg}/50k-3clpcs-0pcs-notypes.bg_genes.tsv"
#   hdest="../../results/eqtl_dynamic/linear_dQTL/${dataset}/${agg}/50k-3clpcs-0pcs-notypes.gsea.tsv"
#   
#   echo $test
#   echo $background
#   echo $hallmark
#   echo $hdest
#   python run_enrichment.py "$test" "$background" "$hallmark" "$hdest"
#   
# done

# for ct in "CM" "CF"
# do
#   test="../../results/eqtl_dynamic/ieQTL/bulk/${ct}/50k-5clpcs-0pcs-notypes.egenes.tsv"
#   background="../../results/eqtl_dynamic/ieQTL/bulk/${ct}/50k-5clpcs-0pcs-notypes.bg_genes.tsv"
#   hdest="../../results/eqtl_dynamic/ieQTL/bulk/${ct}/50k-5clpcs-0pcs-notypes.gsea.tsv"
#   
#   echo $test
#   echo $background
#   echo $hallmark
#   echo $hdest
#   python run_enrichment.py "$test" "$background" "$hallmark" "$hdest"
#   
# done
test="../../results/eqtl_dynamic/linear_dQTL/bulk/day/50k-5clpcs-0pcs-notypes.egenes.tsv"
background="../../results/eqtl_dynamic/linear_dQTL/bulk/day/50k-5clpcs-0pcs-notypes.bg_genes.tsv"
hdest="../../results/eqtl_dynamic/linear_dQTL/bulk/day/50k-5clpcs-0pcs-notypes.gsea.tsv"

echo $test
echo $background
echo $hallmark
echo $hdest
python run_enrichment.py "$test" "$background" "$hallmark" "$hdest"
  
