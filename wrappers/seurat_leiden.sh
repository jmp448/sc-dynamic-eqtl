#!/bin/bash
#SBATCH --partition=gilad
#SBATCH --time=6:00:00
#SBATCH --mem=500G
#SBATCH --output=logs/seurat_leiden.out
#SBATCH --error=logs/seurat_leiden.err
#SBATCH --job-name=leiden

module load R python

Rscript ../code/leiden.R
