#!/bin/bash
#PBS -l select=1:ncpus=16:mem=2gb -l place=pack:excl
# set max execution time
#PBS -l walltime=0:40:00
# set the queue
#PBS -q short_cpuQ

module load mpich-3.2
mpirun.actual -n 16 /home/nicola.muraro/project_hpc/Parallel-implementation-of-Gravitational-Search-Algorithm/src/gsa 30 50 1000 0
mpirun.actual -n 16 /home/nicola.muraro/project_hpc/Parallel-implementation-of-Gravitational-Search-Algorithm/src/gsa 30 100 1000 0
mpirun.actual -n 16 /home/nicola.muraro/project_hpc/Parallel-implementation-of-Gravitational-Search-Algorithm/src/gsa 30 200 1000 0
mpirun.actual -n 16 /home/nicola.muraro/project_hpc/Parallel-implementation-of-Gravitational-Search-Algorithm/src/gsa 30 400 1000 0
mpirun.actual -n 16 /home/nicola.muraro/project_hpc/Parallel-implementation-of-Gravitational-Search-Algorithm/src/gsa 30 800 1000 0
mpirun.actual -n 16 /home/nicola.muraro/project_hpc/Parallel-implementation-of-Gravitational-Search-Algorithm/src/gsa 30 1600 1000 0
mpirun.actual -n 16 /home/nicola.muraro/project_hpc/Parallel-implementation-of-Gravitational-Search-Algorithm/src/gsa 30 3200 1000 0