#!/bin/bash

. /etc/profile

# Configurar OpenMP y otras bibliotecas que usan threads
# usando los valores especificados arriba
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK

module use /opt/ohpc/pub/apps/modules/all/
module load cmake
export CC=$HOME/local/gcc-14.2.0/bin/gcc
export CXX=$HOME/local/gcc-14.2.0/bin/g++
module load GSL
module load OpenBLAS



# Lanzar el programa
srun ${PATH_MACRO}/macro_dr/multi_task/M_scheme_N_tasks.sh


