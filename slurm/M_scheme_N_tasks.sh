#!/bin/bash

. /etc/profile

# Configurar OpenMP y otras bibliotecas que usan threads
# usando los valores especificados arriba
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK

#run cluster specific initialization
source ${PATH_MACRO}/clusters/${CLUSTER}.sh
    
# Lanzar el programa
srun ${PATH_MACRO}/macro_dr/multi_task/M_scheme_N_tasks.sh


