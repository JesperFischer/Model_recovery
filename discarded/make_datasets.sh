#!/bin/bash
#SBATCH --job-name=Hierarchical_Intero_power
#SBATCH --cpus-per-task=15
#SBATCH --mem=2000

# Define your R script and its arguments
R_SCRIPT="datasets_slurm.R"
# Run R script with the specified arguments
Rscript $R_SCRIPT