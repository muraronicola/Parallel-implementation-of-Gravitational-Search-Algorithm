#!/bin/bash
#PBS -l select=1:ncpus=2:mem=2gb -l place=pack:excl
# set max execution time
#PBS -l walltime=6:00:00
# set the queue
#PBS -q short_cpuQ

module load mpich-3.2
mpirun.actual -n 2 /home/nicola.muraro/project_hpc/GSA/gca 2 128000 100 0