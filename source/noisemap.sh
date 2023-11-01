#!/bin/sh

#PBS -N noisemap
#PBS -q Smem
#PBS -l select=1:ncpus=56:ompthreads=56
#PBS -o ./noisemap.out
#PBS -j oe

cd $PBS_O_WORKDIR

for ((i=5;i<1000;i++))
do
    ./noisemap $i
done
