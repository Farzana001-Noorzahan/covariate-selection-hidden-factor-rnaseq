#!/bin/bash
#SBATCH --array=1-16
#SBATCH --job-name=analysis
#SBATCH --output=4-Simulation-sva_%A_%a.out
#SBATCH --error=4-Simulation-sva_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH -p main  # Default and most appropriate partition

# Debugging Information
echo "Running task ID: $SLURM_ARRAY_TASK_ID"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURMD_NODENAME"

# Load R module
module load R/4.4

# Run R script
crun.R Rscript 4-Simulation-sva.R