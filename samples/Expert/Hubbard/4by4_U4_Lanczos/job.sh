#!/bin/sh
#QSUB -queue F4cpu
#QSUB -node  1
#QSUB -mpi   1
#QSUB -omp   24
#QSUB -place pack
#QSUB -over false
#PBS -l walltime=24:00:00
#PBS -N HPhi
cd ${PBS_O_WORKDIR}
 . /etc/profile.d/modules.sh
date
 dplace -x2 ./HPhi -e namelist.def
date
