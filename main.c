#include <stdio.h>
#include "functions.h"
#include "test_functions.h"
#include "utility.h"
#include <stdlib.h>
#include <time.h>

int main(int argc, char *argv[]){
    srand(time(NULL));

    int dim = atoi(argv[1]);
    int pop_size = atoi(argv[2]);
    int n_iter = atoi(argv[3]);

    printf("Configuration: dim = %d, pop_size = %d, n_iter = %d\n", dim, pop_size, n_iter);
    printf("\n");

    float* best_agent = gca(sphere, -100, 100, dim, pop_size, n_iter);

    printf("Best agent: ");
    for (int i = 0; i < dim; i++){
        printf("%f ", best_agent[i]);
    }
    printf("\n");

    printf("Best value: %f\n", sphere(best_agent, dim));

    return 0;
}