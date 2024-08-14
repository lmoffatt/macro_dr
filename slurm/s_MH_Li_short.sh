#!/bin/bash

### Configuración del trabajo

### Nombre de la tarea
#SBATCH --job-name=MH_Li

### Cola a usar (gpu, mono, multi)
#SBATCH --partition=short

### Cantidad de nodos a usar		
#SBATCH --nodes=1

### Cores a utilizar por nodo = procesos por nodo * cores por proceso
#SBATCH --ntasks-per-node=4
### Cores por proceso (para MPI+OpenMP)
#SBATCH --cpus-per-task=16

### Tiempo de ejecucion. Formato dias-horas:minutos.
### short:  <= 1 hora
### multi:  <= 2 días
#SBATCH --time 0-1:00

#---------------------------------------------------

# Script que se ejecuta al arrancar el trabajo

# Cargar el entorno del usuario incluyendo la funcionalidad de modules
# No tocar
. /etc/profile

# Configurar OpenMP y otras bibliotecas que usan threads
# usando los valores especificados arriba
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK

export PATH_MACRO_DR=/home/lmoffatt/macro_dr/v8/

# Cargar los módulos para la tarea
module load amdblis
module load amdlibflame
# Lanzar el programa
srun /home/lmoffatt/macro_dr/macro_dr/multi_task/multi_task_slurm.sh

