#!/bin/bash
#SBATCH --time=10:0
#SBATCH --output=logs/eqtl_stable_mtc_queue.out
#SBATCH --error=logs/eqtl_stable_mtc_queue.err

# rm -r logs/stable-eqtl-mtc/
# mkdir logs/stable-eqtl-mtc/

declare -a CisDists=("50k")
declare -a NumSampPCs=(0)
declare -a NumCLPCs=(5)
declare -a Bins=("bin16")

for cdist in ${CisDists[@]}; do
  for nspc in ${NumSampPCs[@]}; do  
    for ncpc in ${NumCLPCs[@]}; do
      sbatch --mem=15G -J dyn-mtc --error=logs/stable-eqtl-mtc/bulk-$cdist-$nspc-$ncpc-F.err --output=logs/stable-eqtl-mtc/bulk-$cdist-$nspc-$ncpc-F.out eqtl_stable_mtc.sh "bulk" $cdist $nspc $ncpc "F" "day" 10 &
      sbatch --mem=15G -J dyn-mtc --error=logs/stable-eqtl-mtc/bulk7-$cdist-$nspc-$ncpc-F.err --output=logs/stable-eqtl-mtc/bulk7-$cdist-$nspc-$ncpc-F.out eqtl_stable_mtc.sh "bulk7" $cdist $nspc $ncpc "F" "day" 10 &
      sbatch --mem=15G -J dyn-mtc --error=logs/stable-eqtl-mtc/pbd-$cdist-$nspc-$ncpc-F.err --output=logs/stable-eqtl-mtc/pbd-$cdist-$nspc-$ncpc-F.out eqtl_stable_mtc.sh "pseudobulk" $cdist $nspc $ncpc "F" "day" 10 &
      sbatch --mem=15G -J dyn-mtc --error=logs/stable-eqtl-mtc/pbd-$cdist-$nspc-$ncpc-F.err --output=logs/stable-eqtl-mtc/pbd-$cdist-$nspc-$ncpc-F.out eqtl_stable_mtc.sh "pseudobulk-cm" $cdist $nspc $ncpc "F" "day" 10 &
      sbatch --mem=15G -J dyn-mtc --error=logs/stable-eqtl-mtc/pbd-$cdist-$nspc-$ncpc-F.err --output=logs/stable-eqtl-mtc/pbd-$cdist-$nspc-$ncpc-F.out eqtl_stable_mtc.sh "pseudobulk-cf" $cdist $nspc $ncpc "F" "day" 10 &
      for b in ${Bins[@]}; do
        sbatch --mem=15G -J dyn-mtc --error=logs/stable-eqtl-mtc/pb-$b-$cdist-$nspc-$ncpc-F.err --output=logs/stable-eqtl-mtc/pb-$b-$cdist-$nspc-$ncpc-F.out eqtl_stable_mtc.sh "pseudobulk-cm" $cdist $nspc $ncpc "F" $b 10 &
        sbatch --mem=15G -J dyn-mtc --error=logs/stable-eqtl-mtc/pb-$b-$cdist-$nspc-$ncpc-F.err --output=logs/stable-eqtl-mtc/pb-$b-$cdist-$nspc-$ncpc-F.out eqtl_stable_mtc.sh "pseudobulk-cf" $cdist $nspc $ncpc "F" $b 10 &
      done
    done
  done
done

wait
