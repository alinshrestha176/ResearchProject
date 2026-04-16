#!/bin/bash
#SBATCH --job-name=ols_lasso_1000_corr   # Job name
#SBATCH --output=ols_lasso_1000_%j.out   # Standard file output
#SBATCH --error=ols_lasso_1000_%j.err    # Standard file error
#SBATCH --partition=Orion           # Partition name (change if needed)
#SBATCH --nodes=1                   # Number of nodes
#SBATCH --ntasks=1                  # Number of tasks (usually 1 for R scripts)
#SBATCH --cpus-per-task=96          # Number of CPU cores
#SBATCH --mem=192G                   # Memory allocation
#SBATCH --time=40:00:00            # Time limit (hh:mm:ss)

module load R                       # Load R module (modify if needed)

Rscript --vanilla ols_lasso_corrected.R  # Run your R script
