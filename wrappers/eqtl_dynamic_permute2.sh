#!/bin/bash
#SBATCH --time=6:0:0

module load R

Rscript ../code/permutation_dynamic2.R "$@"
