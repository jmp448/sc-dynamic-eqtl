#!/bin/bash
#SBATCH --partition=gilad
#SBATCH --time=10:00:00
#SBATCH --mem=300G
#SBATCH --output=logs/seurat_normalize.out
#SBATCH --error=logs/seurat_normalize.err
#SBATCH --job-name=sctransform

module load R

Rscript ../code/normalize_seurat.R
