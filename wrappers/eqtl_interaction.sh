#!/bin/bash
#SBATCH --time=4:0:0

module load R

Rscript ../code/eqtl_interaction.R "$@"
