#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "utility.h"
#include <stdbool.h>

float **initialize_population(float (*target_function)(float *, int), float lb, float ub, int dim, int pop_size)
{
    float **population = allocate_matrix_float(pop_size, dim);

    for (int i = 0; i < pop_size; i++)
    {
        for (int j = 0; j < dim; j++)
        {
            population[i][j] = random_float(lb, ub);
        }
    }

    return population;
}

float *clip_position_agent(float *agent, float lb, float ub, int dim)
{ // Needed in order to constraint the search space (in the bound of the test function)
    for (int i = 0; i < dim; i++)
    {
        if (agent[i] < lb)
        {
            agent[i] = lb;
        }

        if (agent[i] > ub)
        {
            agent[i] = ub;
        }
    }
    return agent;
}

float get_G(float G0, int t, float n_iter)
{
    return G0 * (1 - t / (n_iter + 1)); // the +1 is needed in order to avoid G0 equal to 0 (at the last iteration)
}

float get_best(float *fitness, int pop_size)
{
    float best = 1e20;
    for (int i = 0; i < pop_size; i++)
    {
        if (fitness[i] < best)
        {
            best = fitness[i];
        }
    }
    return best;
}

float get_worst(float *fitness, int pop_size)
{
    float worst = 0;
    for (int i = 0; i < pop_size; i++)
    {
        if (fitness[i] > worst)
        {
            worst = fitness[i];
        }
    }
    return worst;
}

void sort_agents(float *fitness, float **velocity, float **population, float *M, int pop_size, int dim, int *translation_index)
{
    float tmp_float;
    int tmp_int;
    for (int i = 0; i < pop_size; i++)
    {
        translation_index[i] = i;
    }

    for (int i = 0; i < pop_size; i++)
    {
        for (int j = i + 1; j < pop_size; j++)
        {
            if (fitness[i] > fitness[j])
            {
                tmp_float = fitness[i];
                fitness[i] = fitness[j];
                fitness[j] = tmp_float;

                tmp_int = translation_index[i];
                translation_index[i] = translation_index[j];
                translation_index[j] = tmp_int;

                tmp_float = M[i];
                M[i] = M[j];
                M[j] = tmp_float;

                for (int k = 0; k < dim; k++)
                {
                    tmp_float = velocity[i][k];
                    velocity[i][k] = velocity[j][k];
                    velocity[j][k] = tmp_float;

                    tmp_float = population[i][k];
                    population[i][k] = population[j][k];
                    population[j][k] = tmp_float;
                }
            }
        }
    }
}

float getk_best(int pop_size, int t, float n_iter)
{
    return pop_size * (0.1 + (0.9 * (t / n_iter)));
}

bool check_different_element(float *individual1, float *individual2, int dim)
{
    for (int i = 0; i < dim; i++)
    {
        if (individual1[i] != individual2[i])
        {
            return true;
        }
    }
    return false;
}

float **update_accelerations(float *global_M, float *local_M, float **global_population, float **local_population, float **accelerations, int dim, int pop_size, int k_best, int sub_pop_start_index, int rank, int *translation_index)
{
    float R;
    float **Forces = allocate_matrix_float(pop_size, dim);
    float random;
    int indice;
    printf("Update accelerations\n");
    // printf("\n\nINIZIOOO\n\n");
    /*printf("local_population[0][0]: %f\n", local_population[0][0]);
    printf("local_population[1][0]: %f\n", local_population[1][0]);*/

    for (int i = 0; i < pop_size; i++)
    {
        for (int j = 0; j < k_best; j++)
        {
            if (check_different_element(local_population[i], global_population[j], dim))
            {
                indice = i + sub_pop_start_index;
                printf("\ni: %d; j: %d\n", i, j);
                R = 0;
                for (int d = 0; d < dim; d++)
                {
                    R += (local_population[i][d] - global_population[j][d]) * (local_population[i][d] - global_population[j][d]);
                }
                R = sqrt(R);
                printf("R: %f\n", R);
                printf("M[j]: %f\n", global_M[translation_index[j]]);
                printf("indice %d\n", indice);
                printf("population[i][0] %f\n", local_population[i][0]);
                printf("population[j][0] %f\n", global_population[j][0]);
                printf("translation_index[j] %d\n", translation_index[j]);
                printf("Forces[i][0] %f\n", Forces[i][0]);

                for (int d = 0; d < dim; d++)
                {
                    // random = random_float(0, 1);
                    random = 0.5;
                    // printf("random: %f\n", random);
                    Forces[i][d] = Forces[i][d] + random * global_M[translation_index[j]] * (global_population[j][d] - local_population[i][d]) / (R + 1e-20);
                }
                // printf("Forces[i][0] %f\n", Forces[i][0]);
            }
        }
    }
    printf("Forces[0][0]: %f\n", Forces[0][0]);
    printf("Forces[1][0]: %f\n", Forces[1][0]);

    for (int i = 0; i < pop_size; i++)
    {
        for (int d = 0; d < dim; d++)
        {
            if (local_M[i] > 0)
            { // Mia aggiunta
                accelerations[i][d] = Forces[i][d] / local_M[i];
            }
            else
            {
                accelerations[i][d] = 0;
            }
        }
    }
    /*printf("local_M[0]: %f\n", local_M[0]);
    printf("local_M[1]: %f\n", local_M[1]);
    printf("accelerations[0][0]: %f\n", accelerations[0][0]);
    printf("accelerations[1][0]: %f\n", accelerations[1][0]);*/

    // deallocazione delle matrici
    free(Forces);

    return accelerations;
}

float **update_velocity(float **velocity, float **accelerations, float G, int dim, int pop_size)
{
    float random;
    for (int i = 0; i < pop_size; i++)
    {
        for (int d = 0; d < dim; d++)
        {
            // random = random_float(0, 1);
            random = 0.5;
            velocity[i][d] = random * velocity[i][d] + accelerations[i][d];
        }
    }
    return velocity;
}

float **update_position(float **population, float **velocity, int dim, int pop_size)
{
    for (int i = 0; i < pop_size; i++)
    {
        for (int d = 0; d < dim; d++)
        {
            population[i][d] = population[i][d] + velocity[i][d];
        }
    }
    return population;
}

float **get_local_population(float **local_population, float **global_population, int sub_pop_start_index, int local_pop_size, int global_pop_size, int dim)
{
    for (int i = 0; i < local_pop_size; i++)
    {
        for (int j = 0; j < dim; j++)
        {
            local_population[i][j] = global_population[i + sub_pop_start_index][j];
        }
    }
    return local_population;
}

float *gca(float (*target_function)(float *, int), float lb, float ub, int dim, int global_pop_size, int n_iter, int my_rank, int local_pop_size)
{
    // Returns the best agent found by the algorithm

    // target_function: the function to be optimized (x, dim)
    // lb: lower bound of the search space
    // ub: upper bound of the search space
    // dim: dimensions of the search space
    // pop_size: population size
    // n_iter: number of iterations

    float **global_velocity = allocate_matrix_float(global_pop_size, dim); // Each process calculates for its own subpopulation
    float **local_velocity = allocate_matrix_float(local_pop_size, dim);   // Each process calculates for its own subpopulation
    float **accelerations = allocate_matrix_float(local_pop_size, dim);    // Each process calculates for its own subpopulation
    float *local_fitness = allocate_vector_float(local_pop_size);          // Each process calculates for its own subpopulation
    float *global_fitness = allocate_vector_float(global_pop_size);        // Each process calculates for its own subpopulation
    float *local_M = allocate_vector_float(local_pop_size);                // Each process calculates for its own subpopulation
    float *global_M = allocate_vector_float(global_pop_size);              // Each process calculates for its own subpopulation
    float *m = allocate_vector_float(local_pop_size);
    int *translation_index = allocate_vector_int(global_pop_size); // Each process calculates for its own subpopulation
    float G0 = 100;                                                // G0 = 100
    float G;
    float best, worst;
    float sum_m = 0;
    float local_sum = 0;
    float k_best;
    float **global_population = NULL;
    float **local_population = allocate_matrix_float(local_pop_size, dim);
    int sub_pop_start_index = my_rank * local_pop_size;
    int indice;

    float* best_agents = allocate_vector_float(dim);
    float best_score_best_agent = 1e20;
    /*for (int i = 0; i < global_pop_size; i++)
    {
        translation_index[i] = i;
    }*/

    if (my_rank == 0)
    {
        global_population = initialize_population(target_function, lb, ub, dim, global_pop_size);
        /*printf("Initial population:\n");
        for (int i = 0; i< global_pop_size; i++){
            printf("global_population[%d][0]: %f\n", i, global_population[i][0]);
        }
        printf("-------\n");*/
    }
    else
    {
        global_population = allocate_matrix_float(global_pop_size, dim);
    }
    MPI_Scatter(&(global_population[0][0]), local_pop_size * dim, MPI_FLOAT, &(local_population[0][0]), local_pop_size * dim, MPI_FLOAT, 0, MPI_COMM_WORLD);
    /*printf("Scatter done\n");
    printf("my_rank: %d; local_population[0][0]: %f\n", my_rank, local_population[0][0]);
    printf("my_rank: %d; local_population[1][0]: %f\n", my_rank, local_population[1][0]);*/

    // For additional information
    float *convergence_curve = allocate_vector_float(n_iter);
    float *best_agent = allocate_vector_float(dim);
    float local_best_score = 1e20;

    for (int l = 0; l < n_iter; l++)
    {
        // MPI_Allgather(&(sub_population[0][0]), sub_pop_size * dim, MPI_FLOAT, &(global_population[0][0]), sub_pop_size * dim, MPI_FLOAT, MPI_COMM_WORLD);
        if (my_rank == 0)
        {
            printf("\n\n\nIteration: %d; my_rank: %d\n", l, my_rank);
        }

        for (int i = 0; i < local_pop_size; i++)
        {
            local_population[i] = clip_position_agent(local_population[i], lb, ub, dim);
            local_fitness[i] = target_function(local_population[i], dim);
            if (my_rank == 0)
            {
                printf("my_rank: %d; local_fitness[%d]: %f\n", my_rank, i, local_fitness[i]);
            }

            if (local_fitness[i] < local_best_score)
            {
                local_best_score = local_fitness[i];
            }
        }

        MPI_Allgather(&(local_population[0][0]), local_pop_size * dim, MPI_FLOAT, &(global_population[0][0]), local_pop_size * dim, MPI_FLOAT, MPI_COMM_WORLD);
        MPI_Allgather(local_fitness, local_pop_size, MPI_FLOAT, global_fitness, local_pop_size, MPI_FLOAT, MPI_COMM_WORLD);
        MPI_Allgather(&(local_velocity[0][0]), local_pop_size * dim, MPI_FLOAT, &(global_velocity[0][0]), local_pop_size * dim, MPI_FLOAT, MPI_COMM_WORLD);

        if (my_rank == 0)
        {
            printf("\nnmy_rank: %d; local_population[0][0]: %f\n", my_rank, local_population[0][0]);
            printf("my_rank: %d; local_population[1][0]: %f\n", my_rank, local_population[1][0]);

            printf("my_rank: %d; local_fitness[0]: %f\n", my_rank, local_fitness[0]);
            printf("my_rank: %d; local_fitness[1]: %f\n", my_rank, local_fitness[1]);

            printf("my_rank: %d; local_velocity[0][0]: %f\n", my_rank, local_velocity[0][0]);
            printf("my_rank: %d; local_velocity[1][0]: %f\n", my_rank, local_velocity[1][0]);
        }

        sort_agents(global_fitness, global_velocity, global_population, global_M, global_pop_size, dim, translation_index); // Sort the agents based on their fitness
        k_best = getk_best(global_pop_size, l, n_iter);

        if (my_rank == 0 && best_score_best_agent > global_fitness[0])
        {
            best_score_best_agent = global_fitness[0];
            for (int i = 0; i < dim; i++)
            {
                best_agents[i] = global_population[0][i];
            }
        }

        if (my_rank == 0)
        {
            printf("\nSOORTT:\n");
            printf("my_rank: %d; global_population[0][0]: %f\n", my_rank, global_population[0][0]);
            printf("my_rank: %d; global_population[1][0]: %f\n", my_rank, global_population[1][0]);
            printf("my_rank: %d; global_population[2][0]: %f\n", my_rank, global_population[2][0]);
            printf("my_rank: %d; global_population[3][0]: %f\n", my_rank, global_population[3][0]);

            printf("my_rank: %d; global_fitness[0]: %f\n", my_rank, global_fitness[0]);
            printf("my_rank: %d; global_fitness[1]: %f\n", my_rank, global_fitness[1]);
            printf("my_rank: %d; global_fitness[2]: %f\n", my_rank, global_fitness[2]);
            printf("my_rank: %d; global_fitness[3]: %f\n", my_rank, global_fitness[3]);
        }

        // local_population = get_local_population(local_population, global_population, sub_pop_start_index, local_pop_size, global_pop_size, dim);

        // printf("my_rank: %d; local_population[0][0]: %f\n", my_rank, local_population[0][0]);
        // printf("my_rank: %d; local_population[1][0]: %f\n", my_rank, local_population[1][0]);

        /*for (int i = 0; i < global_pop_size; i++){
            printf("my_rank: %d; global_fitness[%d]: %f\n", my_rank, i, global_fitness[i]);
        }*/

        // Update the G constant
        G = get_G(G0, l, n_iter);
        best = get_best(global_fitness, global_pop_size);
        worst = get_worst(global_fitness, global_pop_size);

        // Update the M and m vectors
        local_sum = 0;
        indice = 0;
        for (int i = 0; i < local_pop_size; i++)
        {
            m[i] = (local_fitness[i] - worst) / (best - worst);
            if (m[i] <= 0)
            { // Mia aggiunta
                m[i] = 0;
            }
            local_sum += m[i];
        }

        MPI_Allreduce(&local_sum, &sum_m, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
        // printf("my_rank: %d; local_population[0][0]: %f\n", my_rank, local_population[0][0]);

        for (int i = 0; i < local_pop_size; i++)
        {
            local_M[i] = m[i] / sum_m;
        }

        /*printf("Before sending local_M: \n");
        printf("%f\n", local_M[0]);
        printf("%f\n", local_M[1]);
        MPI_Allgather(local_M, local_pop_size, MPI_FLOAT, global_M, local_pop_size, MPI_FLOAT, MPI_COMM_WORLD);

        printf("After sending local_M: \n");
        printf("%f\n", local_M[0]);
        printf("%f\n", local_M[1]);*/

        /*printf("\n\n\nmy_rank: %d; local_M[0]: %f\n", my_rank, local_M[0]);
        printf("my_rank: %d; local_M[0]: %f\n", my_rank, local_M[1]);


        printf("\nmy_rank: %d; global_M[0]: %f\n", my_rank, global_M[0]);
        printf("my_rank: %d; global_M[1]: %f\n", my_rank, global_M[1]);
        printf("my_rank: %d; global_M[2]: %f\n", my_rank, global_M[2]);
        printf("my_rank: %d; global_M[3]: %f\n", my_rank, global_M[3]);*/

        MPI_Allgather(local_M, local_pop_size, MPI_FLOAT, global_M, local_pop_size, MPI_FLOAT, MPI_COMM_WORLD);

        if (my_rank == 0)
        {
            printf("\nmy_rank: %d; local_M[0]: %f\n", my_rank, local_M[0]);
            printf("my_rank: %d; local_M[1]: %f\n", my_rank, local_M[1]);

            printf("\nmy_rank: %d; global_M[0]: %f\n", my_rank, global_M[0]);
            printf("my_rank: %d; global_M[1]: %f\n", my_rank, global_M[1]);
            printf("my_rank: %d; global_M[2]: %f\n", my_rank, global_M[2]);
            printf("my_rank: %d; global_M[3]: %f\n", my_rank, global_M[3]);
        }

        accelerations = update_accelerations(global_M, local_M, global_population, local_population, accelerations, dim, local_pop_size, k_best, sub_pop_start_index, my_rank, translation_index);

        if (my_rank == 0)
        {
            printf("\nmy_rank: %d; accelerations[0][0]: %f\n", my_rank, accelerations[0][0]);
            printf("my_rank: %d; accelerations[1][0]: %f\n", my_rank, accelerations[1][0]);
        }

        /*printf("\n\nmy_rank: %d; local_accelerations[0][0]: %f\n", my_rank, accelerations[0][0]);
        printf("my_rank: %d; local_accelerations[1][0]: %f\n", my_rank, accelerations[1][0]);*/

        // printf("my_rank: %d; update_accelerations done\n", my_rank);
        // printf("my_rank: %d; update_accelerations done\n", my_rank);

        local_velocity = update_velocity(local_velocity, accelerations, G, dim, local_pop_size);

        if (my_rank == 0)
        {
            printf("\nmy_rank: %d; local_velocity[0][0]: %f\n", my_rank, local_velocity[0][0]);
            printf("my_rank: %d; local_velocity[1][0]: %f\n", my_rank, local_velocity[1][0]);
        }

        /*printf("\n\nmy_rank: %d; local_velocity[0][0]: %f\n", my_rank, local_velocity[0][0]);
        printf("my_rank: %d; local_velocity[1][0]: %f\n", my_rank, local_velocity[1][0]);*/

        local_population = update_position(local_population, local_velocity, dim, local_pop_size);

        if (my_rank == 0)
        {
            printf("\nmy_rank: %d; local_population[0][0]: %f\n", my_rank, local_population[0][0]);
            printf("my_rank: %d; local_population[1][0]: %f\n", my_rank, local_population[1][0]);
        }

        convergence_curve[l] = local_best_score;

        /*printf("my_rank: %d; local_population[0][0]: %f\n", my_rank, local_population[0][0]);
        printf("my_rank: %d; local_population[1][0]: %f\n", my_rank, local_population[1][0]);*/

        // return local_population[0];

        // printf("Iteration: %d, Best score: %f\n", l, best_score);

        // printf("my_rank: %d; Iteration: %d, Best score: %f\n\n", my_rank, l, local_best_score);
    }

    for (int i = 0; i < local_pop_size; i++)
    {
        local_population[i] = clip_position_agent(local_population[i], lb, ub, dim);
        local_fitness[i] = target_function(local_population[i], dim);

        if (local_fitness[i] < local_best_score)
        {
            local_best_score = local_fitness[i];
            for (int j = 0; j < dim; j++)
            {
                best_agent[j] = local_population[i][j];
            }
        }
    }

    MPI_Allgather(&(local_population[0][0]), local_pop_size * dim, MPI_FLOAT, &(global_population[0][0]), local_pop_size * dim, MPI_FLOAT, MPI_COMM_WORLD);
    MPI_Allgather(local_fitness, local_pop_size, MPI_FLOAT, global_fitness, local_pop_size, MPI_FLOAT, MPI_COMM_WORLD);
    MPI_Allgather(&(local_velocity[0][0]), local_pop_size * dim, MPI_FLOAT, &(global_velocity[0][0]), local_pop_size * dim, MPI_FLOAT, MPI_COMM_WORLD);

    if (my_rank == 0)
    {
        printf("my_rank: %d; local_population[0][0]: %f\n", my_rank, local_population[0][0]);
        printf("my_rank: %d; local_population[0][1]: %f\n", my_rank, local_population[0][1]);
    }

    sort_agents(global_fitness, global_velocity, global_population, global_M, global_pop_size, dim, translation_index); // Sort the agents based on their fitness
    
    if (my_rank == 0 && best_score_best_agent > global_fitness[0])
    {
        best_score_best_agent = global_fitness[0];
        for (int i = 0; i < dim; i++)
        {
            best_agents[i] = global_population[0][i];
        }
    }

    if (my_rank == 0)
    {
        printf("\n---------- DONE ----------\n");
        for (int i = 0; i < global_pop_size; i++)
        {
            printf("my_rank: %d; global_population[%d][0]: %f   global_population[%d][1]: %f \n", my_rank, i, global_population[i][0], i, global_population[i][1]);
        }
        for (int i = 0; i < global_pop_size; i++)
        {
            printf("my_rank: %d; global_fitness[%d]: %f\n", my_rank, i, global_fitness[i]);
        }
    }
    
    //best_agent = global_population[0];
    return best_agents;
}