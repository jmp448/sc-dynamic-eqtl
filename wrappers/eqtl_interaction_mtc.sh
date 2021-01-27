#!/bin/bash
#SBATCH --time=1:0:0

module load R

Rscript ../code/mtc_interaction.R "$@"
