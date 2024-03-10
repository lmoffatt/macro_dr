#!/bin/bash

### Configuración del trabajo

### Nombre de la tarea
#SBATCH --job-name=DR

### Cola a usar (gpu, mono, multi)
#SBATCH --partition=free-rider

### Cantidad de nodos a usar		
#SBATCH --nodes=1

### Cores a utilizar por nodo = procesos por nodo * cores por proceso
#SBATCH --ntasks-per-node=1
### Cores por proceso (para MPI+OpenMP)
#SBATCH --cpus-per-task=64

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
export PATH_MACRO_DR=/nfs/home/lmoffatt/Code/macro_dr/v13/

# Lanzar el programa
srun $PATH_MACRO_DR/macro_dr ../macro_dr/models/scheme_4_inact_NP.txt ../macro_dr/scripts/idealize_experiment.txt "--runIdNamePre= \"v15_DR_\""  "--num_scouts_per_ensemble = get_number(n=32)"  ../macro_dr/scripts/likelihood.txt ../macro_dr/scripts/evidence_thermo_data.txt
