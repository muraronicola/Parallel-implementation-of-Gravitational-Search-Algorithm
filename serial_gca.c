#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "utility.h"


float **serial_initialize_population(float (*target_function)(float*, int), float** velocity, float lb, float ub, int dim, int pop_size, float *fitness, float *M)
{
    float** population = allocate_matrix_float(pop_size, dim);

    for (int i = 0; i < pop_size; i++) {
        for (int j = 0; j < dim; j++) {
            population[i][j] = random_float(lb, ub);
            velocity[i][j] = 0;
        }
        fitness[i] = target_function(population[i], dim);
        M[i] = fitness[i];
    }

    return population;
}

float* serial_clip_position_agent(float* agent, float lb, float ub, int dim) { //Needed in order to constraint the search space (in the bound of the test function)
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


float serial_get_G(float G0, int t, float n_iter){
    return G0 * (1 - t/(n_iter + 1)); //the +1 is needed in order to avoid G0 equal to 0 (at the last iteration)
}

float serial_get_best(float* fitness, int pop_size){
    float best = 1e20;
    for (int i = 0; i < pop_size; i++){
        if (fitness[i] < best){
            best = fitness[i];
        }
    }
    return best;
}

float serial_get_worst(float* fitness, int pop_size){
    float worst = 0;
    for (int i = 0; i < pop_size; i++){
        if (fitness[i] > worst){
            worst = fitness[i];
        }
    }
    return worst;
}

void serial_sort_agents(float* fitness, float** velocity, float** population, float* M, int pop_size, int dim){
    float temp;
    for (int i = 0; i < pop_size; i++){
        for (int j = i + 1; j < pop_size; j++){
            if (fitness[i] > fitness[j]){
                temp = fitness[i];
                fitness[i] = fitness[j];
                fitness[j] = temp;

                temp = M[i];
                M[i] = M[j];
                M[j] = temp;

                for (int k = 0; k < dim; k++){
                    temp = velocity[i][k];
                    velocity[i][k] = velocity[j][k];
                    velocity[j][k] = temp;

                    temp = population[i][k];
                    population[i][k] = population[j][k];
                    population[j][k] = temp;
                }
            }
        }
    }
}


float serial_getk_best(int pop_size, int t, float n_iter){
    return pop_size * (0.1 + (0.9 * (t/n_iter)));
}

float** serial_update_accelearations(float* M, float** population, float** accelerations, int dim, int pop_size, int k_best){
    float R; 
    float **Forces = allocate_matrix_float(pop_size, dim);
    float random;
    for (int i = 0; i < pop_size; i++){
        for (int j = 0; j < k_best; j++){
            if (i != j){
                R = 0;
                for (int d = 0; d < dim; d++){
                    R += (population[i][d] - population[j][d]) * (population[i][d] - population[j][d]);
                }
                R = sqrt(R);

                for (int d = 0; d < dim; d++){
                    random = random_float(0, 1);
                    Forces[i][d] = Forces[i][d] + random * M[j] * (population[j][d] - population[i][d]) / (R + 1e-20);
                }
            }
        }
    }

    for (int i = 0; i < pop_size; i++){
        for (int d = 0; d < dim; d++){
            accelerations[i][d] = Forces[i][d] / M[i];
        }
    }

    //deallocazione delle matrici
    /*for (int i = 0; i < pop_size; i++){ //It was used before
        free(Forces[i]);
    }*/
    free(Forces);

    return accelerations;
}

float ** serial_update_velocity(float** velocity, float** accelerations, float G, int dim, int pop_size){
    float random;
    for (int i = 0; i < pop_size; i++){
        for (int d = 0; d < dim; d++){
            random = random_float(0, 1);
            velocity[i][d] = random * velocity[i][d] + accelerations[i][d];
        }
    }
    return velocity;
}

float** serial_update_position(float** population, float** velocity, int dim, int pop_size){
    for (int i = 0; i < pop_size; i++){
        for (int d = 0; d < dim; d++){
            population[i][d] = population[i][d] + velocity[i][d];
        }
    }
    return population;
}



float* serial_gca(float (*target_function)(float*, int), float lb, float ub, int dim, int pop_size, int n_iter) {
    //Returns the best agent found by the algorithm
    
    // target_function: the function to be optimized (x, dim)
    // lb: lower bound of the search space
    // ub: upper bound of the search space
    // dim: dimensions of the search space
    // pop_size: population size
    // n_iter: number of iterations

    float** velocity = allocate_matrix_float(pop_size, dim);
    float** accelerations = allocate_matrix_float(pop_size, dim);
    float* fitness = allocate_vector_float(pop_size);
    float* M = allocate_vector_float(pop_size);
    float* m = allocate_vector_float(pop_size);
    float G0 = 100; //G0 = 100
    float G;
    float best, worst;
    float sum_m;
    float k_best;

    float** population = serial_initialize_population(target_function, velocity, lb, ub, dim, pop_size, fitness, M);


    //For additional information
    float* convergence_curve = allocate_vector_float(n_iter);
    float* best_agent = allocate_vector_float(dim);
    float best_score = 1e20;



    for (int l = 0; l < n_iter; l++) {
        for (int i = 0; i < pop_size; i++){
            population[i] = serial_clip_position_agent(population[i], lb, ub, dim);
            fitness[i] = target_function(population[i], dim);
            
            if (fitness[i] < best_score) {
                best_score = fitness[i];
                for (int j = 0; j < dim; j++){
                    best_agent[j] = population[i][j];
                }
            }
        }

        serial_sort_agents(fitness, velocity, population, M, pop_size, dim); //Sort the agents based on their fitness
        k_best = serial_getk_best(pop_size, l, n_iter);

        //Update the G constant
        G = serial_get_G(G0, l, n_iter);
        best = serial_get_best(fitness, pop_size);
        worst = serial_get_worst(fitness, pop_size);

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
        accelerations = serial_update_accelearations(M, population, accelerations, dim, pop_size, k_best);
        velocity = serial_update_velocity(velocity, accelerations, G, dim, pop_size);
        population = serial_update_position(population, velocity, dim, pop_size);
    
        convergence_curve[l] = best_score;
    }

    return best_agent;
}