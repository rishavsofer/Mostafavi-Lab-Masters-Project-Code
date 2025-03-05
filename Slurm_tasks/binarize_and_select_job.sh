#!/bin/bash
#SBATCH --job-name=binarize_and_select
#SBATCH --output=binarize_and_select.out
#SBATCH --error=binarize_and_select.err
#SBATCH --time=01:30:00
#SBATCH --mem=40G
#SBATCH --cpus-per-task=8

# Load the Anaconda module
module load anaconda3/cpu/2022.10

# Initialize conda for this shell session
source $(conda info --base)/etc/profile.d/conda.sh

# Activate the conda environment with R installed
conda activate my_R_env

# Run the binarize_and_select.R script (adjust the path if necessary)
Rscript /gpfs/data/mostafavilab/rishav.e.s.dasgupta/R-code/binarize_and_select.R

