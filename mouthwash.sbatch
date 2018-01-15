#!/bin/bash

# This bash script illustrates how to use the sbatch command in SLURM
# to submit a job to the broadwl partition on the RCC cluster. Here is
# an example shell command for running this script using SLURM
# 
#   sbatch mouthwash.sbatch
#

#SBATCH --job-name=mouthwash
#SBATCH --output=mouthwash.out
#SBATCH --error=mouthwash.err
#SBATCH --partition=broadwl
#SBATCH --mem=16G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00

# For reproducibility, add version number to "module load" call.
module load R/3.4.3

# Generate the Makefile targets.
make all