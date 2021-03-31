#!/bin/bash

declare -a CellTypes=("IPSC" "MES" "CMES" "PROG" "CM" "CF")

for ct in ${CellTypes[@]}; do
  cp ../results/eqtl_dynamic/ieQTL/bulk/$ct/50k-5clpcs-0pcs-notypes.mtc.tsv ../tables/ieQTL/${ct}.50k.5clpcs.0pcs.notypes.mtc.tsv
done