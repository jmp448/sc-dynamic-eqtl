#!/bin/bash
#SBATCH --partition=gilad
#SBATCH --time=10:00:00
#SBATCH --mem=300G
#SBATCH --output=logs/seurat_process.out
#SBATCH --error=logs/seurat_process.err
#SBATCH --job-name=umap

module load R

Rscript ../code/process_seurat.R
