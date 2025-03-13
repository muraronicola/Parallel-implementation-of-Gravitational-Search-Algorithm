#!/bin/bash
#PBS -l select=1:ncpus=1:mem=2gb
# set max execution time
#PBS -l walltime=0:30:00
# set the queue
#PBS -q short_cpuQ

module load mpich-3.2
mpirun.actual -n 1 /home/nicola.muraro/project_hpc/GSA/gca 2 32000 100 0


#Reached waaltime dopo 30 minuti, direi che evito questo test