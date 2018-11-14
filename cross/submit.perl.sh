#!/bin/bash

#SBATCH --job-name=HM
#SBATCH --ntasks=1
#SBATCH --partition=gpu

### Environment setup
. /etc/profile

### Environment modules
module load gpu cuda gsl

#OPTIONS += -DCENTROS_MASA_MIN=12.23
#OPTIONS += -DCENTROS_MASA_MAX=13.72
#OPTIONS += -DCENTROS_MASA_MIN=13.72
#OPTIONS += -DCENTROS_MASA_MAX=14.14
#OPTIONS += -DCENTROS_MASA_MIN=14.14
#OPTIONS += -DCENTROS_MASA_MAX=15.34

export LANG=C

### Ejecutar la tarea
srun /mnt/mirta3/marioagustin/halomodelcuda/cross/run_chi.pl 12.23 13.72
#srun /mnt/mirta3/marioagustin/halomodelcuda/cross/run_chi.pl 13.72 14.14
#srun /mnt/mirta3/marioagustin/halomodelcuda/cross/run_chi.pl 14.14 15.34

#srun /mnt/mirta3/marioagustin/halomodelcuda/cross/run.pl 19.0 12.23 13.72
#srun /mnt/mirta3/marioagustin/halomodelcuda/cross/run.pl 19.0 13.72 14.14
#srun /mnt/mirta3/marioagustin/halomodelcuda/cross/run.pl 19.0 14.14 15.34
