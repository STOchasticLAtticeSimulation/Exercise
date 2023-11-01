#!/bin/sh

#PBS -N chaotic
#PBS -q Smem
#PBS -l select=1:ncpus=56:ompthreads=56
#PBS -o ./chaotic.out
#PBS -j oe

cd $PBS_O_WORKDIR

./chaotic 0
