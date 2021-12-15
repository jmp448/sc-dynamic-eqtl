#!/bin/bash
#SBATCH --time=4:0:0
#SBATCH --partition=gilad

module load R

Rscript ../code/mtc_dynamic_permutecells.R "$@"
