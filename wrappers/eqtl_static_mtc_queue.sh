#!/bin/bash
#SBATCH --time=10:0
#SBATCH --mem-per-cpu=5G
#SBATCH --output=logs/eqtl_static_mtc_queue.out
#SBATCH --error=logs/eqtl_static_mtc_queue.err

rm -r logs/static-mtc/
mkdir logs/static-mtc/

declare -a CisDists=("50k")
declare -a NumPCs=(0 1 2 3 4 5)
declare -a BulkDays=("day0" "day1" "day2" "day3" "day4" "day5" "day6" "day7" "day8" "day9" "day10" "day11" "day12" "day13" "day14" "day15")
declare -a PBDays=("day0" "day1" "day3" "day5" "day7" "day11" "day15")
declare -a Types=("iPSC" "meso" "EMT" "cardiomes" "prog" "CM" "EPDC")

for cdist in ${CisDists[@]}; do
  for npc in ${NumPCs[@]}; do
    for d in ${BulkDays[@]}; do
      sbatch -N 1 -n 1 -J bulk-$d-$npc-$cdist --error=logs/static-mtc/bulk-$d-$npc-$cdist.err --output=logs/static-mtc/bulk-$d-$npc-$cdist.out eqtl_static_mtc.sh "bulk" "day" $cdist $npc $d &
    done
    for d in ${PBDays[@]}; do
      sbatch -N 1 -n 1 -J pseudobulk-$d-$npc-$cdist --error=logs/static-mtc/pseudobulk-$d-$npc-$cdist.err --output=logs/static-mtc/pseudobulk-$d-$npc-$cdist.out eqtl_static_mtc.sh "pseudobulk" "day" $cdist $npc $d &
    done
    for ct in ${Types[@]}; do
      sbatch -N 1 -n 1 -J pseudobulk-$ct-$npc-$cdist --error=logs/static-mtc/pseudobulk-$ct-$npc-$cdist.err --output=logs/static-mtc/pseudobulk-$ct-$npc-$cdist.out eqtl_static_mtc.sh "pseudobulk" "type" $cdist $npc $ct &
    done
  done
done

wait