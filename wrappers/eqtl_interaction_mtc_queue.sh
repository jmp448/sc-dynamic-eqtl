#!/bin/bash
#SBATCH --time=10:0
#SBATCH --output=logs/eqtl_interaction_mtc_queue.out
#SBATCH --error=logs/eqtl_interaction_mtc_queue.err

# rm -r logs/interaction-eqtl-mtc/
mkdir logs/interaction-eqtl-mtc/

declare -a Datasets=("bulk")
declare -a CisDists=("50k")
declare -a NumSampPCs=(5 10 20 30)
declare -a NumCLPCs=(5)
declare -a CellTypes=("IPSC" "MES" "CMES" "PROG" "CM" "CF")

for d in ${Datasets[@]}; do
  for cdist in ${CisDists[@]}; do
    for ct in ${CellTypes[@]}; do
      for nspc in ${NumSampPCs[@]}; do
        for ncpc in ${NumCLPCs[@]}; do
          sbatch --mem=15G -J ieqtl-mtc --error=logs/interaction-eqtl-mtc/$d-$cdist-$nspc-$ncpc-F-$ct.err --output=logs/interaction-eqtl-mtc/$d-$cdist-$nspc-$ncpc-F-$ct.out eqtl_interaction_mtc.sh $d $cdist $nspc $ncpc "F" $ct &
        done
      done
    done
  done
done

wait
