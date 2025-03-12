#include <mpi.h>
#include <stdio.h>
#include "parallel_gca.h"
#include "serial_gca.h"
#include "test_functions.h"
#include "utility.h"
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <stdbool.h>

void print_results(double *best_agent, double (*target_function)(double *, int), int dim);

int main(int argc, char *argv[])
{
    MPI_Init(NULL, NULL);

    int comm_sz, my_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    // srand(time(NULL));
    srand(10); // Fix the seed for reproducibility (debug only)

    if (argc != 5)
    {
        printf("Usage: %s <dim> <pop_size> <n_iter> <debug>\n", argv[0]);
        return 1;
    }

    int dim = atoi(argv[1]);
    int pop_size = atoi(argv[2]);
    int n_iter = atoi(argv[3]);
    bool debug = atoi(argv[4]);

    int pop_per_proc = (int)pop_size / comm_sz;
    double *best_agent;
    long seconds, microseconds;
    double elapsed;

    struct timeval begin, end;

    if (comm_sz == 1)
    {
        printf("Configuration: dim = %d, pop_size = %d, n_iter = %d\n", dim, pop_size, n_iter);
        printf("\n");

        gettimeofday(&begin, 0);
        best_agent = serial_gca(sphere, -100, 100, dim, pop_size, n_iter, debug);
        gettimeofday(&end, 0);

        seconds = end.tv_sec - begin.tv_sec;
        microseconds = end.tv_usec - begin.tv_usec;
        elapsed = seconds + microseconds * 1e-6;

        printf("----------------------------------------\n");
        printf("Serial GCA\n");
        printf("----------------------------------------\n");
        print_results(best_agent, sphere, dim);
        printf("Time measured: %.3f seconds.\n", elapsed);
        printf("\n\n----------------------------------------\n");
    }
    else
    {
        srand(10);
        gettimeofday(&begin, 0);
        best_agent = gca(sphere, -100, 100, dim, pop_size, n_iter, my_rank, pop_per_proc, debug);
        gettimeofday(&end, 0);

        if (my_rank == 0)
        {
            printf("Parallel GCA\n");
            printf("----------------------------------------\n");
            seconds = end.tv_sec - begin.tv_sec;
            microseconds = end.tv_usec - begin.tv_usec;
            elapsed = seconds + microseconds * 1e-6;
            print_results(best_agent, sphere, dim);
            printf("Time measured: %.3f seconds.\n", elapsed);
        }
    }

    MPI_Finalize();
    return 0;
}

void print_results(double *best_agent, double (*target_function)(double *, int), int dim)
{
    printf("Best agent: ");
    int i = 0;
    for ( i = 0; i < dim; i++)
    {
        printf("%f ", best_agent[i]);
    }
    printf("\n");

    printf("Best value: %f\n", target_function(best_agent, dim));
}