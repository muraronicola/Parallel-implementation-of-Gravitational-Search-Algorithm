#ifndef _UTILITY_H
#define _UTILITY_H

void check_allocation(void *ptr);
double *allocate_vector_double(int n);
double **allocate_matrix_double(int rows, int columns);
double random_double(int lb, int ub);
int *allocate_vector_int(int n);
double round_to_2_decimals(double number);
double* round_to_2_decimals_vector(double* vector, int dim);
double** round_to_2_decimals_matrix(double** matrix, int rows, int columns);

#endif