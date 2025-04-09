#ifndef __PARALLEL_GSA_H__
#define __PARALLEL_GSA_H__

#include "merge_sort.h"
#include "utility.h"
#include "common.h"
#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <float.h>

double* parallel_gsa(double (*target_function)(double*, int), double lb, double ub, int dim, int global_pop_size, int n_iter, int my_rank, int local_pop_size, bool debug, int n_agents, int* dispacement, int*counts, int* dispacement_matrix, int *count_matrix);
void final_sort(double *source_fitness, double **source_population, double *dest_fitness, double **dest_population, int global_pop_size, int dim, int n_agents, int *dispacement, int *counts);
double **update_accelerations(double *M, double **global_population, double **local_population, double **accelerations, int dim, int pop_size, int k_best, double G);


#endif