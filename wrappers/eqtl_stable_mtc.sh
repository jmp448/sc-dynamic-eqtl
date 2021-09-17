#!/bin/bash
#SBATCH --time=30:0

module load R

Rscript ../code/mtc_stable.R "$@"
