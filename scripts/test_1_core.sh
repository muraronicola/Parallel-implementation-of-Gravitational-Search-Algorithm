#!/bin/bash
#PBS -l select=1:ncpus=1:mem=8gb
# set max execution time
#PBS -l walltime=0:20:00
# set the queue
#PBS -q short_cpuQ

module load mpich-3.2
mpirun.actual -n 1 /home/nicola.muraro/project_hpc/GSA/gca 2 1000 100 0
mpirun.actual -n 1 /home/nicola.muraro/project_hpc/GSA/gca 2 2000 100 0
mpirun.actual -n 1 /home/nicola.muraro/project_hpc/GSA/gca 2 4000 100 0
mpirun.actual -n 1 /home/nicola.muraro/project_hpc/GSA/gca 2 8000 100 0
mpirun.actual -n 1 /home/nicola.muraro/project_hpc/GSA/gca 2 16000 100 0