#!/bin/bash
#SBATCH --job-name=merge_gencode
#SBATCH --output=merge_gencode.out
#SBATCH --error=merge_gencode.err
#SBATCH --time=01:30:00
#SBATCH --mem=20G
#SBATCH --cpus-per-task=8

# Load the Anaconda module
module load anaconda3/cpu/2022.10

# Initialize conda for this shell session
source $(conda info --base)/etc/profile.d/conda.sh

# Activate the conda environment with R installed
conda activate my_R_env

# Run the merge_gencode.R script (update the path if necessary)
Rscript /gpfs/data/mostafavilab/rishav.e.s.dasgupta/R-code/merge_gencode.R

