#!/bin/bash
#SBATCH --partition=gilad
#SBATCH --time=10:00:00
#SBATCH --mem=300G
#SBATCH --output=logs/seurat_process.out
#SBATCH --error=logs/seurat_process.err
#SBATCH --job-name=umap2

module load R

Rscript ../code/cell_type_reprocess.R
