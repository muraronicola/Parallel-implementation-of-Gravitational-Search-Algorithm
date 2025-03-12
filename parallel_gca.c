#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "utility.h"
#include <stdbool.h>
#include <float.h>

double **initialize_population(double (*target_function)(double *, int), double lb, double ub, int dim, int pop_size)
{
    double **population = allocate_matrix_double(pop_size, dim);
    int i = 0;
    for (i = 0; i < pop_size; i++)
    {
        for (int j = 0; j < dim; j++)
        {
            population[i][j] = random_double(lb, ub);
        }
    }

    return population;
}

double *clip_position_agent(double *agent, double lb, double ub, int dim)
{ // Needed in order to constraint the search space (in the bound of the test function)
    int i = 0;
    for (i = 0; i < dim; i++)
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

double get_G(double G0, int t, double n_iter)
{
    return G0 * (1 - t / (n_iter + 1)); // the +1 is needed in order to avoid G0 equal to 0 (at the last iteration)
}

double get_best(double *fitness, int pop_size)
{
    /*
    double best = 1e20;
    for (int i = 0; i < pop_size; i++)
    {
        if (fitness[i] < best)
        {
            best = fitness[i];
        }
    }
    return best;*/
    return fitness[0];
}

double get_worst(double *fitness, int pop_size)
{
    /*double worst = 0;
    for (int i = 0; i < pop_size; i++)
    {
        if (fitness[i] > worst)
        {
            worst = fitness[i];
        }
    }
    return worst;*/
    return fitness[pop_size - 1];
}

void initial_sort(double *fitness, double **population, double *M, int pop_size, int dim, int *translation_index)
{
    double tmp_double;
    int tmp_int;
    int i = 0;
    for (i = 0; i < pop_size; i++)
    {
        translation_index[i] = i;
    }

    int j = 0;
    for (i = 0; i < pop_size; i++)
    {
        for (j = i + 1; j < pop_size; j++)
        {
            if (fitness[i] > fitness[j])
            {
                tmp_double = fitness[i];
                fitness[i] = fitness[j];
                fitness[j] = tmp_double;

                tmp_int = translation_index[i];
                translation_index[i] = translation_index[j];
                translation_index[j] = tmp_int;

                tmp_double = M[i];
                M[i] = M[j];
                M[j] = tmp_double;

                for (int k = 0; k < dim; k++)
                {
                    tmp_double = population[i][k];
                    population[i][k] = population[j][k];
                    population[j][k] = tmp_double;
                }
            }
        }
    }
}

void final_sort(double *source_fitness, double **source_population, double *source_M, int *source_translation_index, double *dest_fitness, double **dest_population, double *dest_M, int *dest_translation_index, int global_pop_size, int local_pop_size, int dim, int n_agents)
{
    int *index_agent = allocate_vector_int(n_agents);

    /*printf("n_agents: %d\n", n_agents);
    printf("local_pop_size: %d\n", local_pop_size);

    for (int i = 0; i < global_pop_size; i++)
    {
        printf("source_fitness[%d]: %f\n", i, source_fitness[i]);
    }*/

    int i = 0;
    for (i = 0; i < n_agents; i++)
    {
        index_agent[i] = i * local_pop_size;
        //printf("index_agent[%d]: %d\n", i, index_agent[i]);
    }


    double lowest_fitness;
    int index_lowest_fitness;
    int j = 0;
    int k = 0;
    for (i = 0; i < global_pop_size; i++)
    {
        lowest_fitness = DBL_MAX;
        for (j = 0; j < n_agents; j++)
        {
            if (index_agent[j] < (j + 1) * local_pop_size)
            {
                if (source_fitness[index_agent[j]] < lowest_fitness)
                {
                    lowest_fitness = source_fitness[index_agent[j]];
                    index_lowest_fitness = j;
                }
            }
        }

        dest_fitness[i] = lowest_fitness;
        dest_M[i] = source_M[index_agent[index_lowest_fitness]];
        dest_translation_index[i] = index_agent[index_lowest_fitness];

        for (k = 0; k < dim; k++)
        {
            dest_population[i][k] = source_population[index_agent[index_lowest_fitness]][k];
        }
        index_agent[index_lowest_fitness]++;
    }

    /*for (int i = 0; i < global_pop_size; i++)
    {
        printf("dest_fitness[%d]: %f\n", i, dest_fitness[i]);
    }

    exit(0);*/
}

double getk_best(int pop_size, int t, double n_iter)
{
    return pop_size * (0.1 + (0.9 * (t / n_iter)));
}

bool check_different_element(double *individual1, double *individual2, int dim)
{
    int i = 0;
    for (i = 0; i < dim; i++)
    {
        if (individual1[i] != individual2[i])
        {
            return true;
        }
    }
    return false;
}

double **update_accelerations(double *global_M, double *local_M, double **global_population, double **local_population, double **accelerations, int dim, int pop_size, int k_best, int sub_pop_start_index, int rank, int *translation_index, bool debug)
{
    double R;
    double **Forces = allocate_matrix_double(pop_size, dim);
    double random;
    int indice;

    if (rank == 0 && debug)
    {
        printf("Update accelerations\n");
    }
    // printf("\n\nINIZIOOO\n\n");
    /*printf("local_population[0][0]: %f\n", local_population[0][0]);
    printf("local_population[1][0]: %f\n", local_population[1][0]);*/

    int i = 0;
    int j = 0;
    int d = 0;
    for (i = 0; i < pop_size; i++)
    {
        for (j = 0; j < k_best; j++)
        {
            if (check_different_element(local_population[i], global_population[j], dim)) //Da cambiare con il check di uguaglianza degli indici
            {
                indice = i + sub_pop_start_index;

                if (rank == 0 && debug)
                {
                    printf("\ni: %d; j: %d\n", i, j);
                }
                R = 0;
                for (d = 0; d < dim; d++)
                {
                    R += (local_population[i][d] - global_population[j][d]) * (local_population[i][d] - global_population[j][d]);
                }
                R = sqrt(R);

                if (rank == 0 && debug)
                {
                    printf("R: %f\n", R);
                    printf("M[j]: %f\n", global_M[translation_index[j]]);
                    printf("indice %d\n", indice);
                    printf("population[i][0] %f\n", local_population[i][0]);
                    printf("population[j][0] %f\n", global_population[j][0]);
                    printf("translation_index[j] %d\n", translation_index[j]);
                    printf("Forces[i][0] %f\n", Forces[i][0]);
                }

                for (d = 0; d < dim; d++)
                {
                    // random = random_double(0, 1);
                    random = 0.5;
                    // printf("random: %f\n", random);
                    Forces[i][d] = Forces[i][d] + random * global_M[translation_index[j]] * (global_population[j][d] - local_population[i][d]) / (R + 1e-20);
                }
                // printf("Forces[i][0] %f\n", Forces[i][0]);
            }
        }
    }

    if (rank == 0 && debug)
    {
        for (i = 0; i < pop_size; i++)
        {
            printf("Forces[%d][0]: %f\n", i, Forces[i][0]);
        }
    }

    for (i = 0; i < pop_size; i++)
    {
        for (d = 0; d < dim; d++)
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

double **update_velocity(double **velocity, double **accelerations, double G, int dim, int pop_size)
{
    double random;
    int i = 0;
    int d = 0;
    for (i = 0; i < pop_size; i++)
    {
        for (d = 0; d < dim; d++)
        {
            // random = random_double(0, 1);
            random = 0.5;
            velocity[i][d] = random * velocity[i][d] + accelerations[i][d];
        }
    }
    return velocity;
}

double **update_position(double **population, double **velocity, int dim, int pop_size)
{
    int i = 0;
    int d = 0;
    for (i = 0; i < pop_size; i++)
    {
        for (d = 0; d < dim; d++)
        {
            population[i][d] = population[i][d] + velocity[i][d];
        }
    }
    return population;
}

double **get_local_population(double **local_population, double **global_population, int sub_pop_start_index, int local_pop_size, int global_pop_size, int dim)
{
    int i = 0;
    int j = 0;
    for (i = 0; i < local_pop_size; i++)
    {
        for (j = 0; j < dim; j++)
        {
            local_population[i][j] = global_population[i + sub_pop_start_index][j];
        }
    }
    return local_population;
}

double *gca(double (*target_function)(double *, int), double lb, double ub, int dim, int global_pop_size, int n_iter, int my_rank, int local_pop_size, bool debug)
{
    // Returns the best agent found by the algorithm

    // target_function: the function to be optimized (x, dim)
    // lb: lower bound of the search space
    // ub: upper bound of the search space
    // dim: dimensions of the search space
    // pop_size: population size
    // n_iter: number of iterations

    // double **global_velocity = allocate_matrix_double(global_pop_size, dim); // Each process calculates for its own subpopulation
    double **local_velocity = allocate_matrix_double(local_pop_size, dim); // Each process calculates for its own subpopulation
    double **accelerations = allocate_matrix_double(local_pop_size, dim);  // Each process calculates for its own subpopulation
    double *local_fitness = allocate_vector_double(local_pop_size);        // Each process calculates for its own subpopulation
    double *global_fitness = allocate_vector_double(global_pop_size);      // Each process calculates for its own subpopulation
    double *local_M = allocate_vector_double(local_pop_size);              // Each process calculates for its own subpopulation
    double *global_M = allocate_vector_double(global_pop_size);            // Each process calculates for its own subpopulation
    double *m = allocate_vector_double(local_pop_size);
    int *translation_index = allocate_vector_int(global_pop_size); // Each process calculates for its own subpopulation
    // int *reverse_translation_index = allocate_vector_int(global_pop_size);
    double G0 = 100; // G0 = 100
    double G;
    double best, worst;
    double sum_m = 0;
    double local_sum = 0;
    double k_best;
    double **global_population = NULL;
    double **local_population = allocate_matrix_double(local_pop_size, dim);
    int sub_pop_start_index = my_rank * local_pop_size;
    int indice;
    int n_agents = global_pop_size / local_pop_size;

    double *unsorted_global_fitness = allocate_vector_double(global_pop_size);
    double *unsorted_global_M = allocate_vector_double(global_pop_size);
    double **unsorted_global_population = allocate_matrix_double(global_pop_size, dim);
    int *unsorted_translation_index = allocate_vector_int(global_pop_size); // Each process calculates for its own subpopulation

    double *best_agents = allocate_vector_double(dim);
    double best_score_best_agent = 1e20;

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
        global_population = allocate_matrix_double(global_pop_size, dim);
    }
    MPI_Scatter(&(global_population[0][0]), local_pop_size * dim, MPI_DOUBLE, &(local_population[0][0]), local_pop_size * dim, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    /*printf("Scatter done\n");
    printf("my_rank: %d; local_population[0][0]: %f\n", my_rank, local_population[0][0]);
    printf("my_rank: %d; local_population[1][0]: %f\n", my_rank, local_population[1][0]);*/

    // For additional information
    double *convergence_curve = allocate_vector_double(n_iter);
    double *best_agent = allocate_vector_double(dim);
    double local_best_score = 1e20;

    int l = 0;
    int i = 0;
    int j = 0;
    for (l = 0; l < n_iter; l++)
    {
        // MPI_Allgather(&(sub_population[0][0]), sub_pop_size * dim, MPI_DOUBLE, &(global_population[0][0]), sub_pop_size * dim, MPI_DOUBLE, MPI_COMM_WORLD);
        if (my_rank == 0 && debug)
        {
            printf("\n\n\nIteration: %d; my_rank: %d\n", l, my_rank);
        }

        for (i = 0; i < local_pop_size; i++)
        {
            local_population[i] = clip_position_agent(local_population[i], lb, ub, dim);
            local_fitness[i] = target_function(local_population[i], dim);
            if (my_rank == 0 && debug)
            {
                printf("my_rank: %d; local_fitness[%d]: %f\n", my_rank, i, local_fitness[i]);
            }

            if (local_fitness[i] < local_best_score)
            {
                local_best_score = local_fitness[i];
            }
        }

        initial_sort(local_fitness, local_population, local_M, local_pop_size, dim, unsorted_translation_index);

        MPI_Allgather(&(local_population[0][0]), local_pop_size * dim, MPI_DOUBLE, &(unsorted_global_population[0][0]), local_pop_size * dim, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Allgather(local_fitness, local_pop_size, MPI_DOUBLE, unsorted_global_fitness, local_pop_size, MPI_DOUBLE, MPI_COMM_WORLD);
        // MPI_Allgather(&(local_velocity[0][0]), local_pop_size * dim, MPI_DOUBLE, &(global_velocity[0][0]), local_pop_size * dim, MPI_DOUBLE, MPI_COMM_WORLD);

        if (my_rank == 0 && debug)
        {
            printf("\n");
            for (i = 0; i < local_pop_size; i++)
            {
                printf("my_rank: %d; local_population[%d][0]: %f\n", my_rank, i, local_population[i][0]);
            }

            printf("\n");
            for (i = 0; i < local_pop_size; i++)
            {
                printf("my_rank: %d; local_fitness[%d]: %f\n", my_rank, i, local_fitness[i]);
            }

            printf("\n");
            for (i = 0; i < local_pop_size; i++)
            {
                printf("my_rank: %d; local_velocity[%d][0]: %f\n", my_rank, i, local_velocity[i][0]);
            }
        }

        final_sort(unsorted_global_fitness, unsorted_global_population, unsorted_global_M, unsorted_translation_index, global_fitness, global_population, global_M, translation_index, global_pop_size, local_pop_size, dim, n_agents);
        // sort_agents(global_fitness, global_population, global_M, global_pop_size, dim, translation_index); // Sort the agents based on their fitness
        /*for(int k = 0; k < global_pop_size; k++){
            reverse_translation_index[translation_index[k]] = k;
        }*/

        k_best = getk_best(global_pop_size, l, n_iter);

        if (my_rank == 0 && best_score_best_agent > global_fitness[0])
        {
            best_score_best_agent = global_fitness[0];
            // printf("Best score: %f\n", best_score_best_agent);
            for (i = 0; i < dim; i++)
            {
                best_agents[i] = global_population[0][i];
            }
        }

        if (my_rank == 0 && debug)
        {
            printf("SORT DONE:\n");
            printf("\n");
            for (i = 0; i < global_pop_size; i++)
            {
                printf("my_rank: %d; global_fitness[%d]: %f\n", my_rank, i, global_fitness[i]);
            }

            printf("\n");
            for (i = 0; i < global_pop_size; i++)
            {
                printf("my_rank: %d; global_population[%d][0]: %f\n", my_rank, i, global_population[i][0]);
            }

            printf("\n");
            for (i = 0; i < global_pop_size; i++)
            {
                printf("my_rank: %d; translation_index[%d]: %d\n", my_rank, i, translation_index[i]);
            }
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
        for (i = 0; i < local_pop_size; i++)
        {
            m[i] = (local_fitness[i] - worst) / (best - worst);
            if (m[i] <= 0)
            { // Mia aggiunta
                m[i] = 0;
            }
            local_sum += m[i];
        }

        MPI_Allreduce(&local_sum, &sum_m, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        // printf("my_rank: %d; local_population[0][0]: %f\n", my_rank, local_population[0][0]);

        if (debug && my_rank == 0)
            printf("it: %d,  sum_m: %.15f\n", l, sum_m);

        for (i = 0; i < local_pop_size; i++)
        {
            local_M[i] = m[i] / sum_m;
        }

        /*printf("Before sending local_M: \n");
        printf("%f\n", local_M[0]);
        printf("%f\n", local_M[1]);
        MPI_Allgather(local_M, local_pop_size, MPI_DOUBLE, global_M, local_pop_size, MPI_DOUBLE, MPI_COMM_WORLD);

        printf("After sending local_M: \n");
        printf("%f\n", local_M[0]);
        printf("%f\n", local_M[1]);*/

        /*printf("\n\n\nmy_rank: %d; local_M[0]: %f\n", my_rank, local_M[0]);
        printf("my_rank: %d; local_M[0]: %f\n", my_rank, local_M[1]);


        printf("\nmy_rank: %d; global_M[0]: %f\n", my_rank, global_M[0]);
        printf("my_rank: %d; global_M[1]: %f\n", my_rank, global_M[1]);
        printf("my_rank: %d; global_M[2]: %f\n", my_rank, global_M[2]);
        printf("my_rank: %d; global_M[3]: %f\n", my_rank, global_M[3]);*/

        // MPI_Allgather(local_M, local_pop_size, MPI_DOUBLE, global_M, local_pop_size, MPI_DOUBLE, MPI_COMM_WORLD); //v1
        MPI_Allgather(local_M, local_pop_size, MPI_DOUBLE, global_M, local_pop_size, MPI_DOUBLE, MPI_COMM_WORLD); // Vorrei poter mandare solo i topk...

        if (my_rank == 0 && debug)
        {
            printf("\n");
            for (i = 0; i < global_pop_size; i++)
            {
                printf("my_rank: %d; global_M[%d]: %f\n", my_rank, i, global_M[i]);
            }

            printf("\n");
            for (i = 0; i < local_pop_size; i++)
            {
                printf("my_rank: %d; local_M[%d]: %f\n", my_rank, i, local_M[i]);
            }
        }

        accelerations = update_accelerations(global_M, local_M, global_population, local_population, accelerations, dim, local_pop_size, k_best, sub_pop_start_index, my_rank, translation_index, debug);

        if (my_rank == 0 && debug)
        {
            printf("\n");
            for (i = 0; i < local_pop_size; i++)
            {
                printf("my_rank: %d; accelerations[%d][0]: %f\n", my_rank, i, accelerations[i][0]);
            }
        }

        /*printf("\n\nmy_rank: %d; local_accelerations[0][0]: %f\n", my_rank, accelerations[0][0]);
        printf("my_rank: %d; local_accelerations[1][0]: %f\n", my_rank, accelerations[1][0]);*/

        // printf("my_rank: %d; update_accelerations done\n", my_rank);
        // printf("my_rank: %d; update_accelerations done\n", my_rank);

        local_velocity = update_velocity(local_velocity, accelerations, G, dim, local_pop_size);

        if (my_rank == 0 && debug)
        {
            printf("\n");
            for (i = 0; i < local_pop_size; i++)
            {
                printf("my_rank: %d; local_velocity[%d][0]: %f\n", my_rank, i, local_velocity[i][0]);
            }
        }

        /*printf("\n\nmy_rank: %d; local_velocity[0][0]: %f\n", my_rank, local_velocity[0][0]);
        printf("my_rank: %d; local_velocity[1][0]: %f\n", my_rank, local_velocity[1][0]);*/

        local_population = update_position(local_population, local_velocity, dim, local_pop_size);

        if (my_rank == 0 && debug)
        {
            printf("\n");
            for (i = 0; i < local_pop_size; i++)
            {
                printf("my_rank: %d; local_population[%d][0]: %f\n", my_rank, i, local_population[i][0]);
            }
        }

        convergence_curve[l] = local_best_score;

        /*printf("my_rank: %d; local_population[0][0]: %f\n", my_rank, local_population[0][0]);
        printf("my_rank: %d; local_population[1][0]: %f\n", my_rank, local_population[1][0]);*/

        // return local_population[0];

        // printf("Iteration: %d, Best score: %f\n", l, best_score);

        // printf("my_rank: %d; Iteration: %d, Best score: %f\n\n", my_rank, l, local_best_score);
    }

    for (i = 0; i < local_pop_size; i++)
    {
        local_population[i] = clip_position_agent(local_population[i], lb, ub, dim);
        local_fitness[i] = target_function(local_population[i], dim);

        if (local_fitness[i] < local_best_score)
        {
            local_best_score = local_fitness[i];
            for (j = 0; j < dim; j++)
            {
                best_agent[j] = local_population[i][j];
            }
        }
    }
    initial_sort(local_fitness, local_population, local_M, local_pop_size, dim, unsorted_translation_index);

    MPI_Allgather(&(local_population[0][0]), local_pop_size * dim, MPI_DOUBLE, &(global_population[0][0]), local_pop_size * dim, MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Allgather(local_fitness, local_pop_size, MPI_DOUBLE, global_fitness, local_pop_size, MPI_DOUBLE, MPI_COMM_WORLD);
    // MPI_Allgather(&(local_velocity[0][0]), local_pop_size * dim, MPI_DOUBLE, &(global_velocity[0][0]), local_pop_size * dim, MPI_DOUBLE, MPI_COMM_WORLD);

    if (my_rank == 0 && debug)
    {
        printf("\n");
        for (i = 0; i < local_pop_size; i++)
        {
            printf("my_rank: %d; local_population[%d][0]: %f\n", my_rank, i, local_population[i][0]);
        }
    }

    final_sort(global_fitness, global_population, global_M, unsorted_translation_index, global_fitness, global_population, global_M, translation_index, global_pop_size, local_pop_size, dim, n_agents);

    if (my_rank == 0 && best_score_best_agent > global_fitness[0])
    {
        best_score_best_agent = global_fitness[0];
        for (i = 0; i < dim; i++)
        {
            best_agents[i] = global_population[0][i];
        }
    }

    if (my_rank == 0 && debug)
    {
        printf("\n---------- DONE ----------\n");
        for (i = 0; i < global_pop_size; i++)
        {
            printf("my_rank: %d; global_population[%d][0]: %f   global_population[%d][1]: %f \n", my_rank, i, global_population[i][0], i, global_population[i][1]);
        }
        for (i = 0; i < global_pop_size; i++)
        {
            printf("my_rank: %d; global_fitness[%d]: %f\n", my_rank, i, global_fitness[i]);
        }
    }

    // best_agent = global_population[0];
    return best_agents;
}