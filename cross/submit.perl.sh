#!/bin/bash

#SBATCH --job-name=HM
#SBATCH --gres=gpu:1
#SBATCH --ntasks=1
#SBATCH --partition=gpu

### Environment setup
. /etc/profile

### Environment modules
module load gcc cuda gsl

### Ejecutar la tarea
srun /home/msgro/halomodelcuda/cross/run_chi.pl
