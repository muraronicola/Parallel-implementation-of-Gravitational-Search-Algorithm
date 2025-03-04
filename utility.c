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

float *allocate_vector_float(int n)
{
    float *ptr = (float *)malloc(sizeof(float) * n);
    check_allocation(ptr);

    return ptr;
}

float **allocate_matrix_float(int rows, int columns)
{
    float *continuos_chunk = (float *)malloc(rows*columns*sizeof(float));
    check_allocation(continuos_chunk);

    float **mat = (float **)malloc(rows * sizeof(float *));
    check_allocation(mat);

    for (int i=0; i<rows; i++)
        mat[i] = &(continuos_chunk[columns*i]);

    return mat;
}

float random_float(int lb, int ub){
    float val = (((float)rand()) /((float)RAND_MAX))*(ub - lb) + lb;
    return val;
}
