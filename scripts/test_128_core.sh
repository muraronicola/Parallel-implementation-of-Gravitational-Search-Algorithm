#!/bin/bash
#PBS -l select=2:ncpus=64:mem=2gb:excl
# set max execution time
#PBS -l walltime=0:10:00
# set the queue
#PBS -q short_cpuQ

module load mpich-3.2
mpirun.actual -n 128 /home/nicola.muraro/project_hpc/GSA/gca 2 1000 100 0
mpirun.actual -n 128 /home/nicola.muraro/project_hpc/GSA/gca 2 2000 100 0
mpirun.actual -n 128 /home/nicola.muraro/project_hpc/GSA/gca 2 4000 100 0
mpirun.actual -n 128 /home/nicola.muraro/project_hpc/GSA/gca 2 8000 100 0
mpirun.actual -n 128 /home/nicola.muraro/project_hpc/GSA/gca 2 16000 100 0
mpirun.actual -n 128 /home/nicola.muraro/project_hpc/GSA/gca 2 32000 100 0
mpirun.actual -n 128 /home/nicola.muraro/project_hpc/GSA/gca 2 64000 100 0
