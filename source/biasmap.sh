#!/bin/sh

#$ -S /bin/sh
#$ -cwd
#$ -V
#$ -q all.q
#$ -N biasmap
#$ -o job_out
#$ -e job_out
#$ -pe OpenMP 36

source /opt/intel/bin/compilervars.sh intel64
export OMP_NUM_THREADS=$NSLOTS
./biasmap
