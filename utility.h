#ifndef _UTILITY_H
#define _UTILITY_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

void check_allocation(void *ptr);
double *allocate_vector_double(int n);
double **allocate_matrix_double(int rows, int columns);
double random_double(int lb, int ub);
int *allocate_vector_int(int n);
double round_to_2_decimals(double number);
double* round_to_2_decimals_vector(double* vector, int dim);
double** round_to_2_decimals_matrix(double** matrix, int rows, int columns);
bool check_different_agents(double *agent_a, double* agent_b, int dim);
void get_displacements_and_counts(int *displacement, int *counts, int *dispacement_matrix, int *count_matrix, int comm_sz, int my_rank, int* pop_per_proc, int remainder, int dim);

#endif