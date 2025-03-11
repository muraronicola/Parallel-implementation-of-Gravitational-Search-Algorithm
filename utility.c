#include <stdio.h>
#include <stdlib.h>

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
    double *continuos_chunk = (double *)malloc(rows*columns*sizeof(double));
    check_allocation(continuos_chunk);
    for (int i=0; i<rows*columns; i++)
        continuos_chunk[i] = 0;

    double **mat = (double **)malloc(rows * sizeof(double *));
    check_allocation(mat);

    for (int i=0; i<rows; i++)
        mat[i] = &(continuos_chunk[columns*i]);

    return mat;
}

double random_double(int lb, int ub){
    double val = (((double)rand()) /((double)RAND_MAX))*(ub - lb) + lb;
    return val;
}
