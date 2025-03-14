#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void check_allocation(void *ptr)
{
    if (ptr == NULL)
    {
        printf("Memory allocation failed\n");
        exit(1);
    }
}

double *allocate_vector_double(int n)
{
    double *ptr = (double *)malloc(sizeof(double) * n);
    check_allocation(ptr);

    return ptr;
}

int *allocate_vector_int(int n)
{
    int *ptr = (int *)malloc(sizeof(int) * n);
    check_allocation(ptr);

    return ptr;
}

double **allocate_matrix_double(int rows, int columns)
{
    double *continuos_chunk = (double *)malloc(rows * columns * sizeof(double));
    check_allocation(continuos_chunk);
    int i = 0;
    for (i = 0; i < rows * columns; i++)
        continuos_chunk[i] = 0;

    double **mat = (double **)malloc(rows * sizeof(double *));
    check_allocation(mat);

    for (i = 0; i < rows; i++)
        mat[i] = &(continuos_chunk[columns * i]);

    return mat;
}

double random_double(int lb, int ub)
{
    double val = (((double)rand()) / ((double)RAND_MAX)) * (ub - lb) + lb;
    return val;
}


double round_to_2_decimals(double number)
{
    return roundf(number * 10000) / 10000;
}

double* round_to_2_decimals_vector(double* vector, int dim)
{
    for (int i = 0; i < dim; i++)
    {
        vector[i] = round_to_2_decimals(vector[i]);
    }
    return vector;
}

double** round_to_2_decimals_matrix(double** matrix, int rows, int columns)
{
    for (int i = 0; i < rows; i++)
    {
       for (int j = 0; j < columns; j++)
       {
            matrix[i][j] = round_to_2_decimals(matrix[i][j]);
       }
    }
    return matrix;
}