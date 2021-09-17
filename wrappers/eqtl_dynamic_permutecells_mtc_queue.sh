#!/bin/bash
#SBATCH --time=10:0
#SBATCH --output=logs/eqtl_dynamic_permutecells_mtc_queue.out
#SBATCH --error=logs/eqtl_dynamic_permutecells_mtc_queue.err

# rm -r logs/dynamic-eqtl-calling/
# mkdir logs/dynamic-eqtl-calling/

declare -a CisDists=("50k")
declare -a NumSampPCs=(0)
declare -a NumCLPCs=(5)
declare -a Bins=("bin16")

for cdist in ${CisDists[@]}; do
  for nspc in ${NumSampPCs[@]}; do
    for ncpc in ${NumCLPCs[@]}; do
      for b in ${Bins[@]}; do
        sbatch --mem=15G -J p_mtc --error=logs/dynamic-eqtl-permutecells/pb-cm-$b-$cdist-$nspc-$ncpc-F-mtc.err --output=logs/dynamic-eqtl-permutecells/pb-cm-$b-$cdist-$nspc-$ncpc-F-mtc.out eqtl_dynamic_permutecells_mtc.sh "pseudobulk-cm" $cdist $nspc $ncpc "F" $b 10 100 &
        # sbatch --mem=15G -J dyn-$c --error=logs/dynamic-eqtl-permutecells/pb-cf-$b-$cdist-$nspc-$ncpc-F-$c-$s.err --output=logs/dynamic-eqtl-permutecells/pb-cf-$b-$cdist-$nspc-$ncpc-F-$c-$s.out eqtl_dynamic_permutecells.sh "pseudobulk-cf" $cdist $nspc $ncpc "F" $b $c 10 $s &
      done
    done
  done
done


wait
