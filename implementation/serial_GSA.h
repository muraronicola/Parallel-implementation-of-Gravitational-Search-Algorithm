#ifndef _SERIAL_GSA_H
#define _SERIAL_GSA_H
#include <stdbool.h>

double* serial_gca(double (*target_function)(double*, int), double lb, double ub, int dim, int pop_size, int n_iter, bool debug);
#endif