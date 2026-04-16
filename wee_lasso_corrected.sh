#!/bin/bash
#SBATCH --job-name=wee_lasso_corr   # Job name
#SBATCH --output=wee_lasso_1000_.out   # Standard file output
#SBATCH --error=wee_lasso_1000_%j.err    # Standard file error
#SBATCH --partition=Orion           # Partition name (change if needed)
#SBATCH --nodes=1                   # Number of nodes
#SBATCH --ntasks=1                  # Number of tasks (usually 1 for R scripts)
#SBATCH --cpus-per-task=96          # Number of CPU cores
#SBATCH --mem=192G                   # Memory allocation
#SBATCH --time=50:00:00            # Time limit (hh:mm:ss)

module load R                       # Load R module (modify if needed)

Rscript --vanilla wee_lasso_diagnosis.R  # Run your R script
