#!/bin/bash
#PBS -l select=1:ncpus=8:mem=2gb
# set max execution time
#PBS -l walltime=0:10:00
# set the queue
#PBS -q short_cpuQ

module load mpich-3.2
mpirun.actual -n 8 /home/nicola.muraro/project_hpc/GSA/gca 2 16000 100 0