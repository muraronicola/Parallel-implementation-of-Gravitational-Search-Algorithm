#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include "parallel_gca.h"
#include "serial_GSA.h"
#include "test_functions.h"
#include "utility.h"
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <stdbool.h>
#include <unistd.h>

void print_results(double *best_agent, double (*target_function)(double *, int), int dim);

int main(int argc, char *argv[])
{
    MPI_Init(NULL, NULL);

    int comm_sz, my_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    if (argc != 5)
    {
        printf("Usage: %s <dim> <pop_size> <n_iter> <debug>\n", argv[0]);
        return 1;
    }

    int dim = atoi(argv[1]);
    int pop_size = atoi(argv[2]);
    int n_iter = atoi(argv[3]);
    bool debug = atoi(argv[4]);

    int pop_per_proc = floor(pop_size / comm_sz);
    int remainder = pop_size % comm_sz;

    unsigned int seed;

    double *best_agent;
    double t1, t2, final_time, best_value;

    if (comm_sz == 1)
    {
        seed = time(NULL);
        srand(seed);

        t1 = MPI_Wtime();
        best_agent = serial_gca(sphere, -1000, 1000, dim, pop_size, n_iter, debug);
        t2 = MPI_Wtime();
        final_time = t2 - t1;

        if (debug)
        {
            printf("Configuration: dim = %d, pop_size = %d, n_iter = %d\n", dim, pop_size, n_iter);
            printf("\n");
            printf("----------------------------------------\n");
            printf("Serial GCA\n");
            printf("----------------------------------------\n");
            print_results(best_agent, sphere, dim);
            printf("Time measured: %.15f seconds.\n", final_time);
            printf("\n\n----------------------------------------\n");
        }
    }
    else
    {
        seed = time(NULL) * (my_rank + 1) * getpid();
        srand(seed);

        t1 = MPI_Wtime();
        int *displacement = (int *)malloc(comm_sz * sizeof(int));
        int *counts = (int *)malloc(comm_sz * sizeof(int));
        int *dispacement_matrix = (int *)malloc(comm_sz * sizeof(int));
        int *count_matrix = (int *)malloc(comm_sz * sizeof(int));
        int i;

        int this_displacement = 0;
        int counter_displacement = 0;

        displacement[0] = 0;
        dispacement_matrix[0] = 0;
        counts[0] = pop_per_proc;
        count_matrix[0] = pop_per_proc;

        if (remainder > 0)
        {
            counts[0]++;
            count_matrix[0]++;
        }
        counter_displacement = counts[0];

        count_matrix[0] *= dim;

        for (i = 1; i < comm_sz; i++)
        {
            this_displacement = pop_per_proc;
            counts[i] = pop_per_proc;
            count_matrix[i] = pop_per_proc;

            if (i < remainder)
            {
                counts[i]++;
                this_displacement++;
                count_matrix[i]++;
            }
            displacement[i] = counter_displacement;
            count_matrix[i] *= dim;
            dispacement_matrix[i] = displacement[i] * dim;
            counter_displacement += counts[i];
        }

        if (my_rank < remainder)
        {
            pop_per_proc++;
        }

        best_agent = gca(sphere, -1000, 1000, dim, pop_size, n_iter, my_rank, pop_per_proc, debug, comm_sz, displacement, counts, dispacement_matrix, count_matrix);
        t2 = MPI_Wtime();
        final_time = t2 - t1;

        if (my_rank == 0)
        {
            if (debug)
            {
                printf("Configuration: dim = %d, pop_size = %d, n_iter = %d\n", dim, pop_size, n_iter);
                printf("\n");
                printf("----------------------------------------\n");
                printf("Parallel GCA\n");
                printf("----------------------------------------\n");
                print_results(best_agent, sphere, dim);
                printf("Time measured: %.15f seconds.\n", final_time);
                printf("\n\n----------------------------------------\n");
            }
        }
    }

    //comm_sz;dim;pop_size;n_iter;seed;time;best_value

    if (my_rank == 0)
    {
        best_value = sphere(best_agent, dim);
        printf("%d;%d;%d;%d;%u;%.15f;%.15f\n", comm_sz, dim, pop_size, n_iter, seed, final_time, best_value);
    }

    MPI_Finalize();
    return 0;
}

void print_results(double *best_agent, double (*target_function)(double *, int), int dim)
{
    printf("Best agent: ");
    int i = 0;
    for (i = 0; i < dim; i++)
    {
        printf("%.15f ", best_agent[i]);
    }
    printf("\n");

    printf("Best value: %.15f\n", target_function(best_agent, dim));
}