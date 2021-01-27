#!/bin/bash
#SBATCH --time=10:0
#SBATCH --mem-per-cpu=15G
#SBATCH --ntasks=32
#SBATCH --output=logs/eqtl_dynamic_permute_queue.out
#SBATCH --error=logs/eqtl_dynamic_permute_queue.err

# rm -r logs/dynamic-eqtl-permute/
mkdir logs/dynamic-eqtl-permute/

declare -a CisDists=("50k")
declare -a NumSampPCs=(0 5 10 20)
declare -a NumCLPCs=(0 3 5 10)

for cdist in ${CisDists[@]}; do
  for nspc in ${NumSampPCs[@]}; do
    for ncpc in ${NumCLPCs[@]}; do
      sbatch -N 1 -n 1 --exclusive -J dyn-b --error=logs/dynamic-eqtl-permute/bulk-$cdist-$nspc-$ncpc.err --output=logs/dynamic-eqtl-permute/bulk-$cdist-$nspc-$ncpc.out eqtl_dynamic_permute.sh "bulk" $cdist $nspc $ncpc "day" &
      #sbatch -N 1 -n 1 --exclusive -J dyn-b7 --error=logs/dynamic-eqtl-calling/bulk7-$cdist-$nspc-$ncpc.err --output=logs/dynamic-eqtl-calling/bulk7-$cdist-$nspc-$ncpc.out eqtl_dynamic_permute.sh "bulk7" $cdist $nspc $ncpc "day" &
      #sbatch -N 1 -n 1 --exclusive -J dyn-pbd --error=logs/dynamic-eqtl-calling/pbd-$cdist-$nspc-$ncpc.err --output=logs/dynamic-eqtl-calling/pbd-$cdist-$nspc-$ncpc.out eqtl_dynamic_permute.sh "pseudobulk" $cdist $nspc $ncpc "day" &
      sbatch -N 1 -n 1 --exclusive -J dyn --error=logs/dynamic-eqtl-permute/pbt-$cdist-$nspc-$ncpc.err --output=logs/dynamic-eqtl-permute/pbt-$cdist-$nspc-$ncpc.out eqtl_dynamic_permute.sh "pseudobulk" $cdist $nspc $ncpc "cmbin" &
      #sbatch -N 1 -n 1 --exclusive -J dyn --error=logs/dynamic-eqtl-calling/pbt-$cdist-$nspc-$ncpc.err --output=logs/dynamic-eqtl-calling/pbt-$cdist-$nspc-$ncpc.out eqtl_dynamic_permute.sh "pseudobulk" $cdist $nspc $ncpc "epdcbin" &
    done
  done
done

wait
