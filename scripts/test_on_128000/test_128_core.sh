#!/bin/bash
#PBS -l select=2:ncpus=64:mem=2gb -l place=excl
# set max execution time
#PBS -l walltime=0:10:00
# set the queue
#PBS -q short_cpuQ

module load mpich-3.2
mpirun.actual -n 128 /home/nicola.muraro/project_hpc/GSA/gca 2 128000 100 0
