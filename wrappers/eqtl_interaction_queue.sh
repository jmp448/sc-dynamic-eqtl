#!/bin/bash
#SBATCH --time=10:0
#SBATCH --output=logs/eqtl_interaction_queue.out
#SBATCH --error=logs/eqtl_interaction_queue.err

# rm -r logs/interaction-eqtl-calling/
mkdir logs/interaction-eqtl-calling/

declare -a Datasets=("bulk")
declare -a CisDists=("50k")
declare -a NumSampPCs=(30)
declare -a NumCLPCs=(5)
declare -a CellTypes=("IPSC" "MES" "CMES" "PROG" "CM" "CF")
declare -a Chunks=(1 2 3 4 5 6 7 8 9 10)

for c in ${Chunks[@]}; do
  for d in ${Datasets[@]}; do
    for cdist in ${CisDists[@]}; do
      for ct in ${CellTypes[@]}; do
        for nspc in ${NumSampPCs[@]}; do
          for ncpc in ${NumCLPCs[@]}; do
            sbatch --mem=15G -J ieqtl-$c --error=logs/interaction-eqtl-calling/$d-$cdist-$nspc-$ncpc-F-$ct-$c.err --output=logs/interaction-eqtl-calling/$d-$cdist-$nspc-$ncpc-F-$ct-$c.out eqtl_interaction.sh $d $cdist $nspc $ncpc "F" $ct $c 10 &
          done
        done
      done
    done
  done
done

wait
