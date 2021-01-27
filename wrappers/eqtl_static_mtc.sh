#!/bin/bash
#SBATCH --time=5:0

module load R

Rscript ../code/mtc_static.R "$@"
