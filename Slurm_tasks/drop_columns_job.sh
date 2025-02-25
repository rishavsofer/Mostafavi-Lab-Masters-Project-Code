#!/bin/bash
#SBATCH --job-name=drop_columns
#SBATCH --output=drop_columns.out
#SBATCH --error=drop_columns.err
#SBATCH --time=01:30:00
#SBATCH --mem=20G
#SBATCH --cpus-per-task=8

# Load the Anaconda module
module load anaconda3/cpu/2022.10

# Initialize conda for this shell session
source $(conda info --base)/etc/profile.d/conda.sh

# Activate the conda environment with R installed
conda activate my_R_env

# Run the drop_columns.R script (update the path if necessary)
Rscript /gpfs/data/mostafavilab/rishav.e.s.dasgupta/R-code/drop_columns.R

