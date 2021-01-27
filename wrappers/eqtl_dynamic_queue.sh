#!/bin/bash
#SBATCH --time=10:0
#SBATCH --mem-per-cpu=15G
#SBATCH --ntasks=80
#SBATCH --output=logs/eqtl_dynamic_queue.out
#SBATCH --error=logs/eqtl_dynamic_queue.err

rm -r logs/dynamic-eqtl-calling/
mkdir logs/dynamic-eqtl-calling/

declare -a CisDists=("50k")
declare -a NumSampPCs=(0 5 10 20)
declare -a NumCLPCs=(0 3 5 10)

for cdist in ${CisDists[@]}; do
  for nspc in ${NumSampPCs[@]}; do
    for ncpc in ${NumCLPCs[@]}; do
      sbatch -N 1 -n 1 --exclusive -J dyn-b --error=logs/dynamic-eqtl-calling/bulk-$cdist-$nspc-$ncpc.err --output=logs/dynamic-eqtl-calling/bulk-$cdist-$nspc-$ncpc.out eqtl_dynamic.sh "bulk" $cdist $nspc $ncpc "day" &
      sbatch -N 1 -n 1 --exclusive -J dyn-b7 --error=logs/dynamic-eqtl-calling/bulk7-$cdist-$nspc-$ncpc.err --output=logs/dynamic-eqtl-calling/bulk7-$cdist-$nspc-$ncpc.out eqtl_dynamic.sh "bulk7" $cdist $nspc $ncpc "day" &
      sbatch -N 1 -n 1 --exclusive -J dyn-pbd --error=logs/dynamic-eqtl-calling/pbd-$cdist-$nspc-$ncpc.err --output=logs/dynamic-eqtl-calling/pbd-$cdist-$nspc-$ncpc.out eqtl_dynamic.sh "pseudobulk" $cdist $nspc $ncpc "day" &
      sbatch -N 1 -n 1 --exclusive -J dyn --error=logs/dynamic-eqtl-calling/pbt-$cdist-$nspc-$ncpc.err --output=logs/dynamic-eqtl-calling/pbt-$cdist-$nspc-$ncpc.out eqtl_dynamic.sh "pseudobulk" $cdist $nspc $ncpc "cmbin" &
      sbatch -N 1 -n 1 --exclusive -J dyn --error=logs/dynamic-eqtl-calling/pbt-$cdist-$nspc-$ncpc.err --output=logs/dynamic-eqtl-calling/pbt-$cdist-$nspc-$ncpc.out eqtl_dynamic.sh "pseudobulk" $cdist $nspc $ncpc "epdcbin" &
    done
  done
done

wait
