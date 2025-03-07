#ifndef _PARALLEL_GCA_H
#define _PARALLEL_GCA_H

#include <stdbool.h>

float* gca(float (*target_function)(float*, int), float lb, float ub, int dim, int global_pop_size, int n_iter, int my_rank, int local_pop_size, bool debug);

#endif