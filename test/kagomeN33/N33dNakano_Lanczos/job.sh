#!/bin/sh
#QSUB -queue B2fat
#QSUB -node  1
#QSUB -mpi   1
#QSUB -omp   40
#QSUB -place pack
#QSUB -over false
#PBS -l walltime=12:00:00
#PBS -N HPhi
cd ${PBS_O_WORKDIR}
 #. /etc/profile.d/modules.sh
#module list > a
#module list
date
mpijob ./HPhi -e namelist.def
date
