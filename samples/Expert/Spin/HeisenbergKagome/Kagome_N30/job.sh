#!/bin/sh
#QSUB -queue F18acc
#QSUB -node  1
#QSUB -mpi   1
#QSUB -omp   24
#QSUB -place pack
#QSUB -over false
#PBS -l walltime=11:00:00
#PBS -N HPhi
cd ${PBS_O_WORKDIR}
 . /etc/profile.d/modules.sh
#module list > a
#module list
date
#export PBS_NUM_PPN=40
#echo "$OMP_NUM_THREADS"
#echo "$PBS_NUM_PPN"
#dplace -x2 ./ED EDnamelist.def
 dplace -x2 ./HPhi -e namelist.def
date
