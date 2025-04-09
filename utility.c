#include "utility.h"

//Check if the memory allocation was successful
void check_allocation(void *ptr)
{
    if (ptr == NULL)
    {
        printf("Memory allocation failed\n");
        exit(1);
    }
}

//Allocate memory for a vector of doubles
double *allocate_vector_double(int n)
{
    double *ptr = (double *)malloc(sizeof(double) * n);
    check_allocation(ptr);

    return ptr;
}

//Allocate memory for a vector of integers
int *allocate_vector_int(int n)
{
    int *ptr = (int *)malloc(sizeof(int) * n);
    check_allocation(ptr);

    return ptr;
}

//Allocate memory for a matrix of floats
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

//Return a random double between lb and ub
double random_double(int lb, int ub)
{
    double val = (((double)rand()) / ((double)RAND_MAX)) * (ub - lb) + lb;
    return val;
}

bool check_different_agents(double *agent_a, double* agent_b, int dim)
{
    int i = 0;
    for (i = 0; i < dim; i++)
    {
        if (agent_a[i] != agent_b[i])
        {
            return true;
        }
    }
    return false;
}

//Gets the number of agents and displacements for each process
void get_displacements_and_counts(int *displacement, int *counts, int *dispacement_matrix, int *count_matrix, int comm_sz, int my_rank, int* pop_per_proc, int remainder, int dim)
{
    int i = 0, counter_displacement = 0;

    displacement[0] = 0;
    dispacement_matrix[0] = 0;
    counts[0] = (*pop_per_proc);
    count_matrix[0] = (*pop_per_proc);

    if (remainder > 0)
    {
        counts[0]++;
        count_matrix[0]++;
    }
    counter_displacement = counts[0];

    count_matrix[0] *= dim;

    for (i = 1; i < comm_sz; i++)
    {
        counts[i] = (*pop_per_proc);
        count_matrix[i] = (*pop_per_proc);

        if (i < remainder)
        {
            counts[i]++;
            count_matrix[i]++;
        }
        displacement[i] = counter_displacement;
        count_matrix[i] *= dim;
        dispacement_matrix[i] = displacement[i] * dim;
        counter_displacement += counts[i];
    }

    if (my_rank < remainder)
    {
        (*pop_per_proc)++;
    }
}