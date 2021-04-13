#!/bin/bash
#SBATCH --time=30:0

module load R

Rscript ../code/permutation_mtc_dynamic3.R "$@"
