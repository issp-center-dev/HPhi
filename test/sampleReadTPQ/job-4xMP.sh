#$ -S /bin/bash
#$ -pe mpich 96
#$ -binding pe linear:96
#$ -cwd
#$ -m n
#$ -v OMP_NUM_THREADS=24
#$ -v P4_RSHCOMMAND=ssh
#$ -v MPICH_PROCESS_GROUP=no
#$ -v CONV_RSH=ssh

echo "module load intel/compiler/64/16.0.1/2016.1.150 intel/mkl/64/11.3.1/2016.1.150 mpich/ge/intel/64/3.2" > ~/.jobenv-$JOB_ID
. ~/.bashrc

printenv

cat $TMPDIR/machines | sed -e 's|:.*||g' > $JOB_NAME.m$JOB_ID

mpirun -machinefile $JOB_NAME.m$JOB_ID -np 4 ./HPhi -e namelistsp.def

/bin/rm -f ~/.jobenv-$JOB_ID
