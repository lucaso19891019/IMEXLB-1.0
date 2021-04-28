#!/bin/bash

#SBATCH --job-name=geng_test
#SBATCH --account=IMEXLBM
#SBATCH --partition=bdwall
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=36
#SBATCH --output=geng_test.out
#SBATCH --error=geng_test.error
#SBATCH --time=01:00:00

# Setup My Environment
module load intel-parallel-studio/cluster.2018.4-ztml34f
export I_MPI_FABRICS=shm:tmi
export OMP_NUM_THREADS=1
export OMP_TARGET_OFFLOAD="DEFAULT"

# Run My Program
srun -n 72 ./run