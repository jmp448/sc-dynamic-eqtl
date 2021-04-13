#!/bin/bash
#SBATCH --partition=gilad
#SBATCH --time=10:00:00
#SBATCH --mem=300G
#SBATCH --output=logs/seurat2h5ad.out
#SBATCH --error=logs/seurat2h5ad.err
#SBATCH --job-name=scloom

module load R/4.0.0
module load hdf5_hl python

Rscript ../code/seurat2h5ad.R
