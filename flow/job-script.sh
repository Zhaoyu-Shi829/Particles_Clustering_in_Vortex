#!/bin/bash
#SBATCH --account=nn9191k
#SBATCH --nodes=9
#SBATCH --job-name=r10f2_22_2
#SBATCH --time=20:00:00
#SBATCH --ntasks-per-node=32

## Recommended safety settings:
set -o errexit # Make bash exit on any error
set -o nounset # Treat unset variables as errors

## Software modules
module reset system
module load intel/2018a
module load HDF5/1.10.1-intel-2018a
module list             # List loaded modules, for easier debugging

# Run simulations
mpirun ./mg41.exe > term_r10f2_22_2.txt 2>&1
mv mglet-perf-report.txt perf-f2_0.txt

