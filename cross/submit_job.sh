#!/bin/bash

#SBATCH --job-name=HM
#SBATCH --gres=gpu:1
#SBATCH --ntasks=1
#SBATCH --partition=gpu
      


### Environment setup
. /etc/profile

### Environment modules

module load compilers/gcc/4.8.2
module load cuda/6.5
module load libs/gsl/1.15-gcc

### Ejecutar la tarea
srun ./hm.cono 0
