#ifndef _UTILITY_H
#define _UTILITY_H

void check_allocation(void *ptr);
double *allocate_vector_double(int n);
double **allocate_matrix_double(int rows, int columns);
double random_double(int lb, int ub);
int *allocate_vector_int(int n);

#endif