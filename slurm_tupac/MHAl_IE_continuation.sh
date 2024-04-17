#!/bin/bash

### Configuración del trabajo

### Nombre de la tarea
#SBATCH --job-name=MHAlIC

### Cola a usar (gpu, mono, multi)
#SBATCH --partition=free-rider

### Cantidad de nodos a usar		
#SBATCH --nodes=1

### Cores a utilizar por nodo = procesos por nodo * cores por proceso
#SBATCH --ntasks-per-node=16
### Cores por proceso (para MPI+OpenMP)
#SBATCH --cpus-per-task=4

### Tiempo de ejecucion. Formato dias-horas:minutos.
### short:  <= 1 hora
### multi:  <= 2 días
#SBATCH --time 2-0:00

#---------------------------------------------------

# Script que se ejecuta al arrancar el trabajo

# Cargar el entorno del usuario incluyendo la funcionalidad de modules
# No tocar
. /etc/profile

# Configurar OpenMP y otras bibliotecas que usan threads
# usando los valores especificados arriba
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK

# configurar el path
export PATH_MACRO_DR=/nfs/home/lmoffatt/Code/macro_dr/v20/

# Lanzar el programa
srun /nfs/home/lmoffatt/Code/macro_dr/macro_dr/multi_task/MH_Al_IE_continuation.sh

