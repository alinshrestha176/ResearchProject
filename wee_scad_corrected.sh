#!/bin/bash
#SBATCH --job-name=wee_scad_1000_500
#SBATCH --output=wee_scad(1000)_%j.out
#SBATCH --error=wee_scad(1000)_%j.err
#SBATCH --partition=Orion
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=96           # 96 cores!
#SBATCH --mem=192G                   # Scale memory with cores
#SBATCH --time=100:00:00              # ~32 hours expected

module load R

Rscript --vanilla wee_scad_corrected.R
