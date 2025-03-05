#include <mpi.h>
#include <stdio.h>
#include "parallel_gca.h"
#include "serial_gca.h"
#include "test_functions.h"
#include "utility.h"
#include <stdlib.h>
#include <time.h>


void print_results(float* best_agent, float (*target_function)(float*, int), int dim);


int main(int argc, char *argv[]){
    MPI_Init(NULL, NULL);

    int comm_sz, my_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    //srand(time(NULL));
    srand(10); //Fix the seed for reproducibility (debug only)

    if (argc != 4){
        printf("Usage: %s <dim> <pop_size> <n_iter>\n", argv[0]);
        return 1;
    }

    int dim = atoi(argv[1]);
    int pop_size = atoi(argv[2]);
    int n_iter = atoi(argv[3]);
    int pop_per_proc = (int) pop_size/comm_sz;
    float* best_agent;


    best_agent = serial_gca(sphere, -100, 100, dim, pop_size, n_iter);
    
    
    if (my_rank == 0){
        printf("Configuration: dim = %d, pop_size = %d, n_iter = %d\n", dim, pop_size, n_iter);
        printf("\n");
    
        printf("----------------------------------------\n");
        printf("Serial GCA\n");
        printf("----------------------------------------\n");
        print_results(best_agent, sphere, dim);
        printf("\n\n----------------------------------------\n");
        printf("Parallel GCA\n");
        printf("----------------------------------------\n");
    }



    best_agent = gca(sphere, -100, 100, dim, pop_size, n_iter, my_rank, pop_per_proc);
    if (my_rank == 0){
        print_results(best_agent, sphere, dim);
    }

    MPI_Finalize();
    return 0;
}


void print_results(float* best_agent, float (*target_function)(float*, int), int dim){
    printf("Best agent: ");
    for (int i = 0; i < dim; i++){
        printf("%f ", best_agent[i]);
    }
    printf("\n");

    printf("Best value: %f\n", target_function(best_agent, dim));
}