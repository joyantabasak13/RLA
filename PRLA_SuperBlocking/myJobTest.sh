#!/bin/bash
#SBATCH -N 5
#SBATCH --partition=hi-core                         # Name of Partition
#SBATCH --ntasks=10                              # Request 256 CPU cores
#SBATCH --time=01:30:00                              # Job should run for up to 1.5 hours (for example)

#SBATCH --partition=epycpriority         # Name of Partition
#SBATCH --constraint='epyc128'
#SBATCH --nodelist=cn485

mpicxx -std=c++17 -O3  -o helloWorld clusterTest.cpp 
srun --mpi=openmpi ./helloWorld    # Replace with your application's commands
# mpirun ./hello
ds7_1M
ds1_50k
ds11_5M_fl
mpicxx -std=c++17 -O3  -o cluster_rla cluster_rla.cpp
srun --mpi=openmpi ./cluster_rla ds7_1M
srun --mpi=openmpi
#SBATCH --nodelist=cn485

#SBATCH -N 6                        # request up to 16 node
#SBATCH --partition=hi-core         # Name of Partition
#SBATCH --ntasks=6                  # Request upto 384 CPU cores
#SBATCH --time=01:30:00             # request upto 6 hours
