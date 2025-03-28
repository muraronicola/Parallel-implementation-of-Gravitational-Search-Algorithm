#!/bin/bash
#PBS -l select=1:ncpus=4:mem=2gb -l place=pack:excl
# set max execution time
#PBS -l walltime=1:00:00
# set the queue
#PBS -q short_cpuQ

module load mpich-3.2
mpirun.actual -n 4 /home/nicola.muraro/project_hpc/GSA/gca 2 1000 100 0
mpirun.actual -n 4 /home/nicola.muraro/project_hpc/GSA/gca 2 2000 100 0
mpirun.actual -n 4 /home/nicola.muraro/project_hpc/GSA/gca 2 4000 100 0
mpirun.actual -n 4 /home/nicola.muraro/project_hpc/GSA/gca 2 8000 100 0
mpirun.actual -n 4 /home/nicola.muraro/project_hpc/GSA/gca 2 16000 100 0
mpirun.actual -n 4 /home/nicola.muraro/project_hpc/GSA/gca 2 32000 100 0
mpirun.actual -n 4 /home/nicola.muraro/project_hpc/GSA/gca 2 64000 100 0