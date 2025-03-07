#ifndef _SERIAL_GCA_H
#define _SERIAL_GCA_H
#include <stdbool.h>

float* serial_gca(float (*target_function)(float*, int), float lb, float ub, int dim, int pop_size, int n_iter, bool debug);
#endif