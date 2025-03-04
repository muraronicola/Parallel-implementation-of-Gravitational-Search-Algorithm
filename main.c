#include <mpi.h>
#include <stdio.h>
#include "functions.h"
#include "test_functions.h"
#include "utility.h"
#include <stdlib.h>
#include <time.h>

int main(int argc, char *argv[]){
    srand(time(NULL));

    if (argc != 4){
        printf("Usage: %s <dim> <pop_size> <n_iter>\n", argv[0]);
        return 1;
    }

    int dim = atoi(argv[1]);
    int pop_size = atoi(argv[2]);
    int n_iter = atoi(argv[3]);

    int comm_sz;
    int my_rank;

    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    printf("Configuration: dim = %d, pop_size = %d, n_iter = %d\n", dim, pop_size, n_iter);
    printf("\n");

    float* best_agent = gca(sphere, -100, 100, dim, pop_size, n_iter);

    printf("Best agent: ");
    for (int i = 0; i < dim; i++){
        printf("%f ", best_agent[i]);
    }
    printf("\n");

    printf("Best value: %f\n", sphere(best_agent, dim));

    MPI_Finalize();
    return 0;
}