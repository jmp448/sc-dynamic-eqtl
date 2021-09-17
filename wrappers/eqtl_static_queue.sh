#!/bin/bash
#SBATCH --time=10:0
#SBATCH --mem-per-cpu=3G
#SBATCH --ntasks=300  
#SBATCH --output=logs/eqtl_static_queue.out
#SBATCH --error=logs/eqtl_static_queue.err

rm -r logs/static-eqtl-calling/
mkdir logs/static-eqtl-calling/

declare -a CisDists=("50k")
declare -a NumPCs=(3 5)
# declare -a NumPCs=(0)
declare -a BulkDays=("day0" "day1" "day2" "day3" "day4" "day5" "day6" "day7" "day8" "day9" "day10" "day11" "day12" "day13" "day14" "day15")
# declare -a BulkDays=("day0")
declare -a PBDays=("day0" "day1" "day3" "day5" "day7" "day11" "day15")
# declare -a PBDays=("day0")
declare -a Types=("IPSC" "MES" "CMES" "PROG" "CM" "CF" "UNK") 
# declare -a Types=("IPSC") 
# declare -a Bins=("bin0" "bin1" "bin2" "bin3" "bin4" "bin5" "bin6" "bin7" "bin8" "bin9" "bin10" "bin11" "bin12" "bin13" "bin14" "bin15")
declare -a Chunks=(1 2 3 4 5)
# declare -a Chunks=(1)

for c in ${Chunks[@]}; do
  for cdist in ${CisDists[@]}; do
    for npc in ${NumPCs[@]}; do
      for d in ${BulkDays[@]}; do
        sbatch -N 1 -n 1 -J bulk-$c --error=logs/static-eqtl-calling/bulk-$d-$npc-$cdist-$c.err --output=logs/static-eqtl-calling/bulk-$d-$npc-$cdist-$c.out eqtl_static.sh "bulk" "day" $cdist $npc $d $c 5 &
      done
      for d in ${PBDays[@]}; do
        sbatch -N 1 -n 1 -J pbd-$c --error=logs/static-eqtl-calling/pseudobulk-$d-$npc-$cdist-$c.err --output=logs/static-eqtl-calling/pseudobulk-$d-$npc-$cdist-$c.out eqtl_static.sh "pseudobulk" "day" $cdist $npc $d $c 5 &
      done
      for ct in ${Types[@]}; do
        sbatch -N 1 -n 1 -J pbt-$c --error=logs/static-eqtl-calling/pseudobulk-$ct-$npc-$cdist-$c.err --output=logs/static-eqtl-calling/pseudobulk-$ct-$npc-$cdist-$c.out eqtl_static.sh "pseudobulk" "type" $cdist $npc $ct $c 5 &
      done
      # for b in ${Bins[@]}; do
      #   sbatch -N 1 -n 1 -J static --error=logs/static-eqtl-calling/pseudobulk-cm-$b-$npc-$cdist.err --output=logs/static-eqtl-calling/pseudobulk-cm-$b-$npc-$cdist.out eqtl_static.sh "pseudobulk-cm" "bin16" $cdist $npc $b &
      #   sbatch -N 1 -n 1 -J static --error=logs/static-eqtl-calling/pseudobulk-cf-$b-$npc-$cdist.err --output=logs/static-eqtl-calling/pseudobulk-cf-$b-$npc-$cdist.out eqtl_static.sh "pseudobulk-cf" "bin16" $cdist $npc $b
      # done
    done
  done
done

wait
