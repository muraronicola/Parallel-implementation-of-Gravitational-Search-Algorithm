#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "utility.h"


float **initialize_population(float (*target_function)(float*, int), float lb, float ub, int dim, int pop_size)
{
    float** population = allocate_matrix_float(pop_size, dim);

    for (int i = 0; i < pop_size; i++) {
        for (int j = 0; j < dim; j++) {
            population[i][j] = random_float(lb, ub);
        }
    }

    return population;
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

void sort_agents(float* fitness, float** velocity, float** population, float* M, int pop_size, int dim){
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


float getk_best(int pop_size, int t, float n_iter){
    return pop_size * (0.1 + (0.9 * (t/n_iter)));
}

float** update_accelerations(float* M, float** population, float** accelerations, int dim, int pop_size, int k_best, int sub_pop_start_index){
    float R; 
    float **Forces = allocate_matrix_float(pop_size, dim);
    float random;
    int indice;
    for (int i = sub_pop_start_index; i < pop_size ; i++){
        for (int j = 0; j < k_best; j++){
            if (i != j){
                indice = j + sub_pop_start_index;
                R = 0;
                for (int d = 0; d < dim; d++){
                    R += (population[i][d] - population[j][d]) * (population[i][d] - population[j][d]);
                }
                R = sqrt(R);

                for (int d = 0; d < dim; d++){
                    random = random_float(0, 1);
                    Forces[i][d] = Forces[i][d] + random * M[indice] * (population[j][d] - population[i][d]) / (R + 1e-20);
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
    free(Forces);

    return accelerations;
}

float ** update_velocity(float** velocity, float** accelerations, float G, int dim, int pop_size){
    float random;
    for (int i = 0; i < pop_size; i++){
        for (int d = 0; d < dim; d++){
            random = random_float(0, 1);
            velocity[i][d] = random * velocity[i][d] + accelerations[i][d];
        }
    }
    return velocity;
}

float** update_position(float** population, float** velocity, int dim, int pop_size){
    for (int i = 0; i < pop_size; i++){
        for (int d = 0; d < dim; d++){
            population[i][d] = population[i][d] + velocity[i][d];
        }
    }
    return population;
}



float* gca(float (*target_function)(float*, int), float lb, float ub, int dim, int global_pop_size, int n_iter, int my_rank, int local_pop_size) {
    // Returns the best agent found by the algorithm
    
    // target_function: the function to be optimized (x, dim)
    // lb: lower bound of the search space
    // ub: upper bound of the search space
    // dim: dimensions of the search space
    // pop_size: population size
    // n_iter: number of iterations

    float** global_velocity = allocate_matrix_float(global_pop_size, dim); //Each process calculates for its own subpopulation
    float** local_velocity = allocate_matrix_float(local_pop_size, dim); //Each process calculates for its own subpopulation
    float** accelerations = allocate_matrix_float(local_pop_size, dim); //Each process calculates for its own subpopulation
    float* local_fitness = allocate_vector_float(local_pop_size); //Each process calculates for its own subpopulation
    float* global_fitness = allocate_vector_float(global_pop_size); //Each process calculates for its own subpopulation
    float* local_M = allocate_vector_float(local_pop_size); //Each process calculates for its own subpopulation
    float* global_M = allocate_vector_float(global_pop_size); //Each process calculates for its own subpopulation
    float* m = allocate_vector_float(local_pop_size); //Each process calculates for its own subpopulation
    float G0 = 100; //G0 = 100
    float G;
    float best, worst;
    float sum_m = 0;
    float local_sum = 0;
    float k_best;
    float** global_population = NULL;
    float** local_population = allocate_matrix_float(local_pop_size, dim);
    int sub_pop_start_index = my_rank * local_pop_size;
    int indice;

    if (my_rank == 0){
        global_population = initialize_population(target_function, lb, ub, dim, global_pop_size);
    }else{
        global_population = allocate_matrix_float(global_pop_size, dim);
    }
    MPI_Scatter(&(global_population[0][0]), local_pop_size * dim, MPI_FLOAT, &(local_population[0][0]), local_pop_size * dim, MPI_FLOAT, 0, MPI_COMM_WORLD);

    //For additional information
    float* convergence_curve = allocate_vector_float(n_iter);
    float* best_agent = allocate_vector_float(dim);
    float local_best_score = 1e20;

    for (int l = 0; l < n_iter; l++) {
        //MPI_Allgather(&(sub_population[0][0]), sub_pop_size * dim, MPI_FLOAT, &(global_population[0][0]), sub_pop_size * dim, MPI_FLOAT, MPI_COMM_WORLD);

        for (int i = 0; i < local_pop_size; i++){
            local_population[i] = clip_position_agent(local_population[i], lb, ub, dim);
            local_fitness[i] = target_function(local_population[i], dim);
            
            if (local_fitness[i] < local_best_score) {
                local_best_score = local_fitness[i];
                for (int j = 0; j < dim; j++){
                    best_agent[j] = local_population[i][j];
                }
            }
        }

        MPI_Allgather(&(local_population[0][0]), local_pop_size * dim, MPI_FLOAT, &(global_population[0][0]), local_pop_size * dim, MPI_FLOAT, MPI_COMM_WORLD);
        MPI_Allgather(local_fitness, local_pop_size * dim, MPI_FLOAT, global_fitness, local_pop_size * dim, MPI_FLOAT, MPI_COMM_WORLD);
        MPI_Allgather(&(local_velocity[0][0]), local_pop_size * dim, MPI_FLOAT, &(global_velocity[0][0]), local_pop_size * dim, MPI_FLOAT, MPI_COMM_WORLD);
        printf("my_rank: %d; MPI_Allgather done\n", my_rank);

        sort_agents(global_fitness, global_velocity, global_population, global_M, global_pop_size, dim); //Sort the agents based on their fitness
        k_best = getk_best(global_pop_size, l, n_iter);

        //Update the G constant
        G = get_G(G0, l, n_iter);
        best = get_best(global_fitness, global_pop_size);
        worst = get_worst(global_fitness, global_pop_size);

        //Update the M and m vectors
        local_sum = 0;
        indice = 0;
        for (int i = 0; i < local_pop_size; i++){
            indice = i + sub_pop_start_index;
            m[i] = (global_fitness[indice] - worst) / (best - worst);
            local_sum += m[i];
        }

        MPI_Allreduce(&local_sum, &sum_m, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);

        for (int i = 0; i < local_pop_size; i++){
            local_M[i] = m[i] / sum_m;
        }

        MPI_Allgather(local_M, local_pop_size * dim, MPI_FLOAT, global_M, local_pop_size * dim, MPI_FLOAT, MPI_COMM_WORLD);

        //Update the velocity
        //printf("my_rank: %d;   pop[0][0]: %f\n", my_rank, population[0][0]);
        //printf("my_rank: %d;   sub_population[0][0]: %f\n", my_rank, sub_population[0][0]);
        //printf("my_rank: %d;   update_accelerations\n", my_rank);
        //printf("my_rank: %d;   pop[0][0]: %f\n", my_rank, population[0][0]);
        //printf("my_rank: %d;   sub_population[0][0]: %f\n", my_rank, sub_population[0][0]);
        //printf("my_rank: %d;   accelerations[0][0]: %f\n", my_rank, accelerations[0][0]);
        //printf("my_rank: %d;   M[0]: %f\n", my_rank, M[0]);
        //printf("my_rank: %d;   k_best: %f\n", my_rank, k_best);

        accelerations = update_accelerations(global_M, global_population, accelerations, dim, local_pop_size, k_best, sub_pop_start_index);
        
        //printf("my_rank: %d; update_accelerations done\n", my_rank);
        //printf("my_rank: %d; update_accelerations done\n", my_rank);

        local_velocity = update_velocity(local_velocity, accelerations, G, dim, local_pop_size);
        //printf("my_rank: %d; update_velocity done\n", my_rank);
        local_population = update_position(local_population, local_velocity, dim, local_pop_size);
        //printf("my_rank: %d; update_position done\n", my_rank);
    
        convergence_curve[l] = local_best_score;
        //printf("Iteration: %d, Best score: %f\n", l, best_score);

        printf("my_rank: %d; Iteration: %d, Best score: %f\n\n", my_rank, l, local_best_score);
    }

    return best_agent;
}