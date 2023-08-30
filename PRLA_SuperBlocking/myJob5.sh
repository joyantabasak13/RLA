#!/bin/bash
#SBATCH --partition=hi-core         # partition for MPI
#SBATCH -N 1                        # Can request up to 16 node
#SBATCH --ntasks=5                  # Can request upto 384 CPU cores
#SBATCH --time=01:30:00             # Can request upto 6 hours
#SBATCH --nodelist=cn485
#SBATCH --mem=128G

srun --mpi=openmpi ./prla_superBlocking_perfectLoadBalancing ds11_5M_fl
