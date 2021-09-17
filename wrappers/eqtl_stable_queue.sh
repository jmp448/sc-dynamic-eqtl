#!/bin/bash
#SBATCH --time=10:0
#SBATCH --output=logs/eqtl_dynamic_queue.out
#SBATCH --error=logs/eqtl_dynamic_queue.err

# rm -r logs/dynamic-eqtl-calling/
# mkdir logs/dynamic-eqtl-calling/

declare -a CisDists=("50k")
declare -a NumSampPCs=(0)
declare -a NumCLPCs=(6 7 8 9 10)
declare -a Bins=("bin16")
declare -a Chunks=(1 2 3 4 5 6 7 8 9 10)

# without regressing out cell type proportions
for c in ${Chunks[@]}; do
  for cdist in ${CisDists[@]}; do
    for nspc in ${NumSampPCs[@]}; do
      for ncpc in ${NumCLPCs[@]}; do
        sbatch --mem=15G -J dyn-$c --error=logs/dynamic-eqtl-calling/bulk-$cdist-$nspc-$ncpc-F-$c.err --output=logs/dynamic-eqtl-calling/bulk-$cdist-$nspc-$ncpc-F-$c.out eqtl_dynamic.sh "bulk" $cdist $nspc $ncpc "F" "day" $c 10 &
#         sbatch --mem=15G -J dyn-$c --error=logs/dynamic-eqtl-calling/bulk7-$cdist-$nspc-$ncpc-F-$c.err --output=logs/dynamic-eqtl-calling/bulk7-$cdist-$nspc-$ncpc-F-$c.out eqtl_dynamic.sh "bulk7" $cdist $nspc $ncpc "F" "day" $c 10 &
#         sbatch --mem=15G -J dyn-$c --error=logs/dynamic-eqtl-calling/pb-d-$cdist-$nspc-$ncpc-F-$c.err --output=logs/dynamic-eqtl-calling/pb-d-$cdist-$nspc-$ncpc-F-$c.out eqtl_dynamic.sh "pseudobulk" $cdist $nspc $ncpc "F" "day" $c 10 &
#         sbatch --mem=15G -J dyn-$c --error=logs/dynamic-eqtl-calling/pb-cm-d-$cdist-$nspc-$ncpc-F-$c.err --output=logs/dynamic-eqtl-calling/pb-cm-d-$cdist-$nspc-$ncpc-F-$c.out eqtl_dynamic.sh "pseudobulk-cm" $cdist $nspc $ncpc "F" "day" $c 10 &
#         sbatch --mem=15G -J dyn-$c --error=logs/dynamic-eqtl-calling/pb-cf-d-$cdist-$nspc-$ncpc-F-$c.err --output=logs/dynamic-eqtl-calling/pb-cf-d-$cdist-$nspc-$ncpc-F-$c.out eqtl_dynamic.sh "pseudobulk-cf" $cdist $nspc $ncpc "F" "day" $c 10 &
        for b in ${Bins[@]}; do
          sbatch --mem=15G -J dyn-$c --error=logs/dynamic-eqtl-calling/pb-cm-$b-$cdist-$nspc-$ncpc-F-$b-$c.err --output=logs/dynamic-eqtl-calling/pb-cm-$b-$cdist-$nspc-$ncpc-F-$c.out eqtl_dynamic.sh "pseudobulk-cm" $cdist $nspc $ncpc "F" $b $c 10 &
          sbatch --mem=15G -J dyn-$c --error=logs/dynamic-eqtl-calling/pb-cf-$b-$cdist-$nspc-$ncpc-F-$b-$c.err --output=logs/dynamic-eqtl-calling/pb-cf-$b-$cdist-$nspc-$ncpc-F-$c.out eqtl_dynamic.sh "pseudobulk-cf" $cdist $nspc $ncpc "F" $b $c 10 &
        done
      done
    done
  done
done

# with cell type proportion regression
# for c in ${Chunks[@]}; do
#   for cdist in ${CisDists[@]}; do
#     for nspc in ${NumSampPCs[@]}; do
#       for ncpc in ${NumCLPCs[@]}; do
#         # sbatch --mem=15G -J dyn-$c --error=logs/dynamic-eqtl-calling/bulk-$cdist-$nspc-$ncpc-T-$c.err --output=logs/dynamic-eqtl-calling/bulk-$cdist-$nspc-$ncpc-T-$c.out eqtl_dynamic.sh "bulk" $cdist $nspc $ncpc "T" "day" $c 10 &
#         sbatch --mem=15G -J dyn-$c --error=logs/dynamic-eqtl-calling/bulk7-$cdist-$nspc-$ncpc-T-$c.err --output=logs/dynamic-eqtl-calling/bulk7-$cdist-$nspc-$ncpc-T-$c.out eqtl_dynamic.sh "bulk7" $cdist $nspc $ncpc "T" "day" $c 10 &
#         # sbatch --mem=15G -J dyn-$c --error=logs/dynamic-eqtl-calling/pb-d-$cdist-$nspc-$ncpc-T-$c.err --output=logs/dynamic-eqtl-calling/pb-d-$cdist-$nspc-$ncpc-T-$c.out eqtl_dynamic.sh "pseudobulk" $cdist $nspc $ncpc "T" "day" $c 10 &
#       done
#     done
#   done
# done

wait
