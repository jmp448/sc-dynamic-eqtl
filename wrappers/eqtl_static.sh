#!/bin/bash
#SBATCH --time=3:0:0

module load R

Rscript ../code/eqtl_static.R "$@"
