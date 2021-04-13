#!/bin/bash
#SBATCH --partition=gilad
#SBATCH --time=6:00:00
#SBATCH --mem=150G
#SBATCH --output=logs/ti.out
#SBATCH --error=logs/ti.err
#SBATCH --job-name=ti

module load python
source activate singlecell

python ../code/ti.py
