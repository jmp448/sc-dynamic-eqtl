#!/bin/bash
#SBATCH --time=2:0:0

module load R/4.0.0

Rscript ../code/mashr_static.R "$@"
