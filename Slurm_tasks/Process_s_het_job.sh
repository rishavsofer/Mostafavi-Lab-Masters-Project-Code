#!/bin/bash
#SBATCH --job-name=process_s_het
#SBATCH --output=process_s_het.out
#SBATCH --error=process_s_het.err
#SBATCH --time=01:30:00
#SBATCH --mem=40G
#SBATCH --cpus-per-task=8

# Load the Anaconda module
module load anaconda3/cpu/2022.10

# Initialize conda for this shell session
source $(conda info --base)/etc/profile.d/conda.sh

# Activate the conda environment with R installed
conda activate my_R_env

# Run the Process_s_het.R script (update the path if necessary)
Rscript /gpfs/data/mostafavilab/rishav.e.s.dasgupta/R-code/Process_s_het.R

