#include <stdio.h>
#include "utility.h"


float **initialize_population(float (*target_function)(float*, int), float** velocity, float lb, float ub, int dim, int pop_size, float *fitness, float *M, float* gBestScore)
{
    float** population = allocate_matrix_float(pop_size, dim);

    for (int i = 0; i < pop_size; i++) {
        for (int j = 0; j < dim; j++) {
            population[i][j] = random_float(lb, ub);
            velocity[i][j] = 0;
        }
        fitness[i] = target_function(population[i], dim);
        M[i] = fitness[i];
        if (fitness[i] < *gBestScore) {
            *gBestScore = fitness[i];
        }
    }
}

float* clip_position_agent(float* agent, float lb, float ub, int dim) { //Needed in order to constraint the search space (in the bound of the test function)
    for (int i = 0; i < dim; i++) {
        if (agent[i] < lb) {
            agent[i] = lb;
        }
        else if (agent[i] > ub) {
            agent[i] = ub;
        }
    }
    return agent;
}


float get_G(float G0, int t, float n_iter){
    return G0 * (1 - t/(n_iter + 1)); //the +1 is needed in order to avoid G0 equal to 0 (at the last iteration)
}

float get_best(float* fitness, int pop_size){
    float best = 1e20;
    for (int i = 0; i < pop_size; i++){
        if (fitness[i] < best){
            best = fitness[i];
        }
    }
    return best;
}

float get_worst(float* fitness, int pop_size){
    float worst = 0;
    for (int i = 0; i < pop_size; i++){
        if (fitness[i] > worst){
            worst = fitness[i];
        }
    }
    return worst;
}




void gca(float (*target_function)(float*, int), float lb, float ub, int dim, int pop_size, int n_iter) {
    // target_function: the function to be optimized (x, dim)
    // lb: lower bound of the search space
    // ub: upper bound of the search space
    // dim: dimensions of the search space
    // pop_size: population size
    // n_iter: number of iterations

    float** velocity = allocate_matrix_float(pop_size, dim);
    float* fitness = allocate_vector_float(pop_size);
    float* M = allocate_vector_float(pop_size);
    float* m = allocate_vector_float(pop_size);
    float gBestScore = 1e20;
    float G0 = 100; //G0 = 100
    float G;
    float best, worst;
    float sum_m;

    float** population = initialize_population(target_function, velocity, lb, ub, dim, pop_size, fitness, M, &gBestScore);


    for (int l = 0; l < n_iter; l++) {
        for (int i = 0; i < pop_size; i++){
            population[i] = clip_position_agent(population[i], lb, ub, dim);
            fitness[i] = target_function(population[i], dim);
        }

        //Update the G constant
        G = get_G(G0, l, n_iter);
        best = get_best(fitness, pop_size);
        worst = get_worst(fitness, pop_size);




        //Update the M and m vectors
        sum_m = 0;
        for (int i = 0; i < pop_size; i++){
            m[i] = (fitness[i] - worst) / (best - worst);
            sum_m += m[i];
        }
        for (int i = 0; i < pop_size; i++){
            M[i] = m[i] / sum_m;
        }
        
        //Update the velocity
        
    }
}