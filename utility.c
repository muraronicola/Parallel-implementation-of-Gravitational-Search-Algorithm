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
    float **mat = (float **)malloc(rows * sizeof(float *));
    check_allocation(mat);

    for (int i = 0; i < rows; i++)
    {
        mat[i] = (float *)malloc(columns * sizeof(float));
        check_allocation(mat[i]);
    }

    return mat;
}


float random_float(int lb, int ub){
    return rand() % (ub + 1 - lb) + lb;
}