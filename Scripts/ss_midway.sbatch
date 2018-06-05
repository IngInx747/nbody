#!/bin/bash
#SBATCH --time=00:10:00
#SBATCH --partition=sandyb
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH --exclusive

module load openmpi

## strong scaling against nodes

# nodes=8 rank-per-node=16 thread-per-rank=1
mpirun -n 128 ./nbody-parallel.exe --threads=1 --particles=102400 --iterations=10

# nodes=8 rank-per-node=1 thread-per-rank=16
mpirun -n 8 ./nbody-parallel.exe --threads=16 --particles=102400 --iterations=10