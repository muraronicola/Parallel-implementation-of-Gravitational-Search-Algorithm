#ifndef _UTILITY_H
#define _UTILITY_H

void check_allocation(void *ptr);
float *allocate_vector_float(int n);
float **allocate_matrix_float(int rows, int columns);
float random_float(int lb, int ub);

#endif