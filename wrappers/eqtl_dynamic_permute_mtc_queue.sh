#!/bin/bash
#SBATCH --time=10:0
#SBATCH --mem-per-cpu=5G
#SBATCH --ntasks=8
#SBATCH --output=logs/eqtl_dynamic_permute_mtc_queue.out
#SBATCH --error=logs/eqtl_dynamic_permute_mtc_queue.err

rm -r logs/dynamic-eqtl-permute-mtc/
mkdir logs/dynamic-eqtl-permute-mtc/

declare -a CisDists=("50k")
declare -a NumSampPCs=(0)
declare -a NumCLPCs=(5)
declare -a Bins=("bin16")

for cdist in ${CisDists[@]}; do
  for nspc in ${NumSampPCs[@]}; do  
    for ncpc in ${NumCLPCs[@]}; do
      sbatch --mem=15G -J dynp-mtc --error=logs/dynamic-eqtl-permute-mtc/bulk-$cdist-$nspc-$ncpc-F.err --output=logs/dynamic-eqtl-permute-mtc/bulk-$cdist-$nspc-$ncpc-F.out eqtl_dynamic_permute_mtc.sh "bulk" $cdist $nspc $ncpc "F" "day" 10 &
      # sbatch --mem=15G -J dynp-mtc --error=logs/dynamic-eqtl-permute-mtc/bulk7-$cdist-$nspc-$ncpc-F.err --output=logs/dynamic-eqtl-permute-mtc/bulk7-$cdist-$nspc-$ncpc-F.out eqtl_dynamic_permute_mtc.sh "bulk7" $cdist $nspc $ncpc "F" "day" 10 &
      # sbatch --mem=15G -J dynp-mtc --error=logs/dynamic-eqtl-permute-mtc/pbd-$cdist-$nspc-$ncpc-F.err --output=logs/dynamic-eqtl-permute-mtc/pbd-$cdist-$nspc-$ncpc-F.out eqtl_dynamic_permute_mtc.sh "pseudobulk" $cdist $nspc $ncpc "F" "day" 10 &
      # sbatch --mem=15G -J dynp-mtc --error=logs/dynamic-eqtl-permute-mtc/pbcmd-$cdist-$nspc-$ncpc-F.err --output=logs/dynamic-eqtl-permute-mtc/pbcmd-$cdist-$nspc-$ncpc-F.out eqtl_dynamic_permute_mtc.sh "pseudobulk-cm" $cdist $nspc $ncpc "F" "day" 10 &
      # sbatch --mem=15G -J dynp-mtc --error=logs/dynamic-eqtl-permute-mtc/pbcfd-$cdist-$nspc-$ncpc-F.err --output=logs/dynamic-eqtl-permute-mtc/pbcfd-$cdist-$nspc-$ncpc-F.out eqtl_dynamic_permute_mtc.sh "pseudobulk-cf" $cdist $nspc $ncpc "F" "day" 10 &
      for b in ${Bins[@]}; do
        sbatch --mem=15G -J dynp-mtc --error=logs/dynamic-eqtl-permute-mtc/pbcm-$b-$cdist-$nspc-$ncpc-F.err --output=logs/dynamic-eqtl-permute-mtc/pbcm-$b-$cdist-$nspc-$ncpc-F.out eqtl_dynamic_permute_mtc.sh "pseudobulk-cm" $cdist $nspc $ncpc "F" $b 10 &
        sbatch --mem=15G -J dynp-mtc --error=logs/dynamic-eqtl-permute-mtc/pbcf-$b-$cdist-$nspc-$ncpc-F.err --output=logs/dynamic-eqtl-permute-mtc/pbcf-$b-$cdist-$nspc-$ncpc-F.out eqtl_dynamic_permute_mtc.sh "pseudobulk-cf" $cdist $nspc $ncpc "F" $b 10 &
      done
    done
  done
done

wait