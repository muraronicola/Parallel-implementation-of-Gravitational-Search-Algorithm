#!/bin/bash
#PBS -l select=1:ncpus=32:mem=2gb -l place=pack:excl
# set max execution time
#PBS -l walltime=0:10:00
# set the queue
#PBS -q short_cpuQ

module load mpich-3.2
mpirun.actual -n 32 /home/nicola.muraro/project_hpc/GSA/gca 30 100 1000 0
mpirun.actual -n 32 /home/nicola.muraro/project_hpc/GSA/gca 30 200 1000 0
mpirun.actual -n 32 /home/nicola.muraro/project_hpc/GSA/gca 30 400 1000 0
mpirun.actual -n 32 /home/nicola.muraro/project_hpc/GSA/gca 30 800 1000 0
mpirun.actual -n 32 /home/nicola.muraro/project_hpc/GSA/gca 30 1600 1000 0
mpirun.actual -n 32 /home/nicola.muraro/project_hpc/GSA/gca 30 3200 1000 0
mpirun.actual -n 32 /home/nicola.muraro/project_hpc/GSA/gca 30 6400 1000 0