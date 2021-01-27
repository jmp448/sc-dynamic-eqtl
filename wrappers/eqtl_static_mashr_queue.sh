#!/bin/bash
#SBATCH --time=10:0
#SBATCH --mem-per-cpu=50G
#SBATCH --ntasks=15
#SBATCH --output=logs/eqtl_static_mashr_queue.out
#SBATCH --error=logs/eqtl_static_mashr_queue.err

rm -r logs/static-mashr/
mkdir logs/static-mashr/

declare -a CisDists=("50k")
declare -a NumPCs=(0 1 2 3 4 5)

for cdist in ${CisDists[@]}; do
  for npc in ${NumPCs[@]}; do
    sbatch -N 1 -n 1 --exclusive -J mashr --error=logs/static-mashr/bulk-$d-$npc-$cdist.err --output=logs/static-mashr/bulk-$d-$npc-$cdist.out eqtl_static_mashr.sh "bulk" "day" $cdist $npc &
    sbatch -N 1 -n 1 --exclusive -J mashr --error=logs/static-mashr/pseudobulk-$d-$npc-$cdist.err --output=logs/static-mashr/pseudobulk-$d-$npc-$cdist.out eqtl_static_mashr.sh "pseudobulk" "day" $cdist $npc &
    sbatch -N 1 -n 1 --exclusive -J mashr --error=logs/static-mashr/pseudobulk-$ct-$npc-$cdist.err --output=logs/static-mashr/pseudobulk-$ct-$npc-$cdist.out eqtl_static_mashr.sh "pseudobulk" "type" $cdist $npc &
  done
done

wait
