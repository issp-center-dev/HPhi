#!/bin/sh
#QSUB -queue i18cpu
#QSUB -node  16
#QSUB -mpi   16
#QSUB -omp   24
#QSUB -place pack
#QSUB -over false
#PBS -l walltime=00:30:00
#PBS -N HPhi
cd ${PBS_O_WORKDIR}
 #. /etc/profile.d/modules.sh
#module list > a
#module list
date
 mpijob ./HPhi -e namelist.def
date
