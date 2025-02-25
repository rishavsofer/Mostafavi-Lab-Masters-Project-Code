#!/bin/bash
#SBATCH --job-name=check_gene_features
#SBATCH --output=check_gene_features_%j.out
#SBATCH --error=check_gene_features_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=00:30:00
#SBATCH --mem=20G
#SBATCH --partition=cpu_short

# Load the Anaconda module
module load anaconda3/cpu/2022.10
# Initialize conda for this shell session
source $(conda info --base)/etc/profile.d/conda.sh
# Activate the conda environment with R installed
conda activate my_R_env

# Define the path to your R script (.rtf file)
R_SCRIPT="/gpfs/data/mostafavilab/rishav.e.s.dasgupta/R-code/checking_gene_features.R"

# Run the R script using Rscript
Rscript "${R_SCRIPT}"

