#ifndef __SERIAL_GSA_H__
#define __SERIAL_GSA_H__

#include "utility.h"
#include "merge_sort.h"
#include "common.h"
#include <stdbool.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

double* serial_gsa(double (*target_function)(double*, int), double lb, double ub, int dim, int pop_size, int n_iter);
double **serial_update_accelearations(double *M, double **population, double **accelerations, int dim, int pop_size, int k_best, double G);

#endif