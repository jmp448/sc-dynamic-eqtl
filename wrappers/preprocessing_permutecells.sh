#!/bin/bash
#SBATCH --time=10:0:0
#SBATCH --mem=100G
#SBATCH --partition=bigmem2
#SBATCH --error=logs/preprocessing_permutecells.err
#SBATCH --output=logs/preprocessing_permutecells.out

module load R

Rscript ../code/preprocessing_permutecells.R
