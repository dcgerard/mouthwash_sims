#!/bin/bash

# This bash script illustrates how to use the sbatch command in SLURM
# to submit a job to the broadwl partition on the RCC cluster. Here is
# an example shell command for running this script using SLURM
# 
#   sbatch mouthwash.sbatch
#
# Note that the "time" option may have to be adjusted; 36 hours may
# not be long enough to complete all the computation.

#SBATCH --job-name=mouthwash
#SBATCH --output=mouthwash.out
#SBATCH --error=mouthwash.err
#SBATCH --partition=broadwl
#SBATCH --mem=16G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=36:00:00

# For reproducibility, add version number to "module load" call.
module load R/3.4.3

# The packages installed here are sometimes not detected
# automatically, so I add this location manually to the package search
# path.
export R_LIBS_USER=~/R/x86_64-redhat-linux-gnu-library/3.4

# Set the number of threads used for matrix operations. This should
# speed up the runtime considerably because there are lots of
# intensive matrix computations.
export OPENBLAS_NUM_THREADS=8

# Generate the Makefile targets.
make all
