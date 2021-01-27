#!/bin/bash
#SBATCH --time=10:0
#SBATCH --mem-per-cpu=5G
#SBATCH --ntasks=8
#SBATCH --output=logs/eqtl_dynamic_mtc_queue.out
#SBATCH --error=logs/eqtl_dynamic_mtc_queue.err

# rm -r logs/dynamic-eqtl-mtc/
# mkdir logs/dynamic-eqtl-mtc/

declare -a CisDists=("50k")
# declare -a NumSampPCs=(0 5 10 20)
declare -a NumSampPCs=(0 5 10 20)
declare -a NumCLPCs=(0 3 5 10)

for cdist in ${CisDists[@]}; do
  for nspc in ${NumSampPCs[@]}; do
    for ncpc in ${NumCLPCs[@]}; do
      # sbatch -N 1 -n 1 --exclusive -J dyn-mtc --error=logs/dynamic-eqtl-mtc/bulk-$cdist-$nspc-$ncpc.err --output=logs/dynamic-eqtl-mtc/bulk-$cdist-$nspc-$ncpc.out eqtl_dynamic_mtc.sh "bulk" $cdist $nspc $ncpc "day" &
      # sbatch -N 1 -n 1 --exclusive -J dyn-mtc --error=logs/dynamic-eqtl-mtc/bulk7-$cdist-$nspc-$ncpc.err --output=logs/dynamic-eqtl-mtc/bulk7-$cdist-$nspc-$ncpc.out eqtl_dynamic_mtc.sh "bulk7" $cdist $nspc $ncpc "day" &
      sbatch -N 1 -n 1 --exclusive -J dyn-mtc --error=logs/dynamic-eqtl-mtc/pbd-$cdist-$nspc-$ncpc.err --output=logs/dynamic-eqtl-mtc/pbd-$cdist-$nspc-$ncpc.out eqtl_dynamic_mtc.sh "pseudobulk" $cdist $nspc $ncpc "epdcbin" &
      sbatch -N 1 -n 1 --exclusive -J dyn-mtc --error=logs/dynamic-eqtl-mtc/pbt-$cdist-$nspc-$ncpc.err --output=logs/dynamic-eqtl-mtc/pbt-$cdist-$nspc-$ncpc.out eqtl_dynamic_mtc.sh "pseudobulk" $cdist $nspc $ncpc "cmbin" &
    done
  done
done

wait
