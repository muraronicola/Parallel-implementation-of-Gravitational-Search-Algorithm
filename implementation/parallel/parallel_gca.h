#ifndef _PARALLEL_GCA_H
#define _PARALLEL_GCA_H

#include <stdbool.h>

double* gca(double (*target_function)(double*, int), double lb, double ub, int dim, int global_pop_size, int n_iter, int my_rank, int local_pop_size, bool debug, int n_agents, int* dispacement, int*counts, int* dispacement_matrix, int *count_matrix);

#endif