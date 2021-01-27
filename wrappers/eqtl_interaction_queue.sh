#!/bin/bash
#SBATCH --time=10:0
#SBATCH --mem-per-cpu=20G
#SBATCH --ntasks=6
#SBATCH --output=logs/eqtl_interaction_queue.out
#SBATCH --error=logs/eqtl_interaction_queue.err

# rm -r logs/interaction-eqtl-calling/
# mkdir logs/interaction-eqtl-calling/

# declare -a Datasets=("bulk" "pseudobulk")
declare -a Datasets=("bulk")
declare -a CisDists=("50k")
# declare -a NumSampPCs=(0 5 10 20)
declare -a NumSampPCs=(5)
# declare -a NumCLPCs=(0 3 5 10)
declare -a NumCLPCs=(5)
declare -a CellTypes=("iPSC" "meso" "cardiomes" "prog" "EPDC" "CM")

for d in ${Datasets[@]}; do
  for cdist in ${CisDists[@]}; do
    for ct in ${CellTypes[@]}; do
      for nspc in ${NumSampPCs[@]}; do
        for ncpc in ${NumCLPCs[@]}; do
          sbatch -N 1 -n 1 --exclusive -J ieqtl --error=logs/interaction-eqtl-calling/$d-$cdist-$nspc-$ncpc-$ct.err --output=logs/interaction-eqtl-calling/$d-$cdist-$nspc-$ncpc-$ct.out eqtl_interaction.sh $d $cdist $nspc $ncpc $ct &
        done
      done
    done
  done
done

wait
