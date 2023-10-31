#!/bin/sh

#PBS -N biasmap
#PBS -q Smem
#PBS -l select=1:ncpus=56:ompthreads=56
#PBS -o ./biasmap.out
#PBS -j oe

cd $PBS_O_WORKDIR

./biasmap
