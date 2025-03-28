#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "utility.h"
#include <stdbool.h>
#include <float.h>
#include "merge_sort.h"

void sort_agents(double *fitness, double **population, double *M, int pop_size, int dim, int *translation_index)
{
    double tmp_double;
    int tmp_int;
    int i = 0;
    int j = 0;
    int k = 0;
    for (i = 0; i < pop_size; i++)
    {
        translation_index[i] = i;
    }

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

                for (k = 0; k < dim; k++)
                {
                    tmp_double = population[i][k];
                    population[i][k] = population[j][k];
                    population[j][k] = tmp_double;
                }
            }
        }
    }
}

double **initialize_population(double (*target_function)(double *, int), double lb, double ub, int dim, int pop_size)
{
    double **population = allocate_matrix_double(pop_size, dim);
    int i = 0;
    int j = 0;
    for (i = 0; i < pop_size; i++)
    {
        for (j = 0; j < dim; j++)
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

double get_G(double G0, int t, int n_iter)
{
    double t_double = (double)t;
    double n_iter_double = (double)n_iter;
    double result = G0 * exp(-20 * (t_double / n_iter_double));
    //printf("G: %f\n", result);
    return result;
    // return G0 * (1 - t / (n_iter + 1)); // the +1 is needed in order to avoid G0 equal to 0 (at the last iteration)
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

void initial_sort(double *fitness, double **population, double *local_fitness_sorted, double **local_population_sorted, int pop_size, int dim, int *local_translation_index)
{
    double tmp_double;
    int tmp_int;
    int i = 0;
    int j = 0;
    for (i = 0; i < pop_size; i++)
    {
        local_translation_index[i] = i;
        // printf("local_translation_index[%d]: %d\n", i, local_translation_index[i]);
    }

    for (i = 0; i < pop_size; i++)
    {
        local_fitness_sorted[i] = fitness[i];
        for (j = 0; j < dim; j++)
        {
            local_population_sorted[i][j] = population[i][j];
        }
    }

    // printf("-----------------\n");
    double temp;
    int k = 0;
    int this_element = 0;
    int next_element = 0;
    for (i = 0; i < pop_size - 1; i++)
    {
        for (j = 0; j < pop_size - i - 1; j++)
        {
            this_element = j;
            next_element = j + 1;

            if (local_fitness_sorted[this_element] > local_fitness_sorted[next_element])
            {
                temp = local_fitness_sorted[this_element];
                local_fitness_sorted[this_element] = local_fitness_sorted[next_element];
                local_fitness_sorted[next_element] = temp;

                /*temp = M[this_element];
                M[this_element] = M[next_element];
                M[next_element] = temp;*/

                tmp_int = local_translation_index[this_element];
                local_translation_index[this_element] = local_translation_index[next_element];
                local_translation_index[next_element] = tmp_int;

                for (k = 0; k < dim; k++)
                {
                    temp = local_population_sorted[this_element][k];
                    local_population_sorted[this_element][k] = local_population_sorted[next_element][k];
                    local_population_sorted[next_element][k] = temp;

                    /*temp = local_velocity[this_element][k];
                    local_velocity[this_element][k] = local_velocity[next_element][k];
                    local_velocity[next_element][k] = temp;*/
                }
            }
        }
    }
}

void final_sort(double *source_fitness, double **source_population, int *unsorted_translation_index, double *dest_fitness, double **dest_population, int *dest_translation_index, int global_pop_size, int local_pop_size, int dim, int n_agents, int *dispacement, int *counts, int *dispacement_matrix, int *count_matrix)
{
    int *index_agent = allocate_vector_int(n_agents);

    /*printf("local_translation_index[0]: %d\n", local_translation_index[0]);
    printf("local_translation_index[1]: %d\n", local_translation_index[1]);
    printf("local_translation_index[2]: %d\n", local_translation_index[2]);
    printf("local_translation_index[3]: %d\n", local_translation_index[3]);
    printf("local_translation_index[4]: %d\n", local_translation_index[4]);



    printf("dest_translation_index[0]: %d\n", dest_translation_index[0]);
    printf("dest_translation_index[1]: %d\n", dest_translation_index[1]);
    printf("dest_translation_index[2]: %d\n", dest_translation_index[2]);
    printf("dest_translation_index[3]: %d\n", dest_translation_index[3]);
    printf("dest_translation_index[4]: %d\n", dest_translation_index[4]);
    exit(0);*/

    /*printf("n_agents: %d\n", n_agents);
    printf("local_pop_size: %d\n", local_pop_size);

    for (int i = 0; i < global_pop_size; i++)
    {
        printf("source_fitness[%d]: %f\n", i, source_fitness[i]);
    }*/

    int i = 0;
    int j = 0;
    for (i = 0; i < n_agents; i++)
    {
        index_agent[i] = dispacement[i];
        // printf("index_agent[%d]: %d\n", i, index_agent[i]);
    }

    /*printf("Initial unsorted_translation_index\n");

    for (i = 0; i < global_pop_size; i++)
    {
        printf("unsorted_translation_index[%d]: %d\n", i, unsorted_translation_index[i]);
    }*/

    double multiplier = 0;
    j = 0;
    int counter_pop = 0;
    for (i = 0; i < global_pop_size; i++)
    {

        /*if (i % local_pop_size == 0 && i != 0)
        {
            multiplier++;
        }
        unsorted_translation_index[i] = unsorted_translation_index[i] + (multiplier * local_pop_size);
        */
        unsorted_translation_index[i] = unsorted_translation_index[i] + dispacement[j];
        counter_pop++;
        if (counter_pop == counts[j])
        {
            j++;
            counter_pop = 0;
        }
    }

    /*printf("unsorted_translation_index\n");

    for (i = 0; i < global_pop_size; i++)
    {
        printf("unsorted_translation_index[%d]: %d\n", i, unsorted_translation_index[i]);
    }*/

    double lowest_fitness;
    int index_lowest_fitness;
    int k = 0;
    for (i = 0; i < global_pop_size; i++)
    {
        lowest_fitness = DBL_MAX;
        for (j = 0; j < n_agents; j++)
        {
            if (index_agent[j] < counts[j] + dispacement[j])
            {
                if (source_fitness[index_agent[j]] < lowest_fitness)
                {
                    lowest_fitness = source_fitness[index_agent[j]];
                    index_lowest_fitness = j;
                }
            }
        }

        dest_fitness[i] = lowest_fitness;
        dest_translation_index[i] = unsorted_translation_index[index_agent[index_lowest_fitness]];

        for (k = 0; k < dim; k++)
        {
            dest_population[i][k] = source_population[index_agent[index_lowest_fitness]][k];
        }
        index_agent[index_lowest_fitness]++;
    }

    /*printf("final dest_translation_index\n");

    for (i = 0; i < global_pop_size; i++)
    {
        printf("dest_translation_index[%d]: %d\n", i, dest_translation_index[i]);
    }

    exit(0);*/
    /*for (int i = 0; i < global_pop_size; i++)
    {
        printf("dest_fitness[%d]: %f\n", i, dest_fitness[i]);
    }

    exit(0);*/
}

double getk_best(int pop_size, int t, int n_iter)
{
    // return pop_size * (0.1 + (0.9 * (t / n_iter)));
    double t_double = (double)t;
    double n_iter_double = (double)n_iter;
    double result = pop_size * ((n_iter_double - t_double) / n_iter_double);
    result = ceil(result);
    //printf("result: %f\n", result);
    return result;
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

double **update_accelerations(double *global_M, double *local_M, double **global_population, double **local_population, double **accelerations, int dim, int pop_size, int k_best, int sub_pop_start_index, int rank, int *translation_index, double G, bool debug)
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
            indice = i + sub_pop_start_index;

            if (indice != translation_index[j]) // Da cambiare con il check di uguaglianza degli indici
            {
                if (rank == 0 && debug)
                {
                    printf("\ni: %d; j: %d\n", i, j);
                }
                R = 0;
                for (d = 0; d < dim; d++)
                {
                    R += pow(local_population[i][d] - global_population[j][d], 2);
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
                    //random = 0.5;
                    // printf("random: %f\n", random);
                    //Forces[i][d] = Forces[i][d] + random * (G * ((global_M[translation_index[j]] * local_M[i]) / (R + 1e-20)) * (global_population[j][d] - local_population[i][d]));
                    Forces[i][d] = Forces[i][d] + (((global_M[translation_index[j]]) / (R + 1e-20)) * (global_population[j][d] - local_population[i][d]));
                    
                }
                // printf("Forces[i][0] %f\n", Forces[i][0]);
            }
        }

        for (d = 0; d < dim; d++)
        {
            //random = random_double(0, 1);
            random = 1;
            Forces[i][d] = random * Forces[i][d] * G;
        }
    }

    if (rank == 0 && debug)
    {
        for (i = 0; i < pop_size; i++)
        {
            printf("my_rank: %d; Forces[%d][0]: %f\n", rank, i, Forces[i][0]);
        }
    }

    for (i = 0; i < pop_size; i++)
    {
        for (d = 0; d < dim; d++)
        {
            accelerations[i][d] = Forces[i][d];
            /*if (Forces[i][d] > 0)
            {
                accelerations[i][d] = Forces[i][d] / local_M[i];
            }
            else
            {
                accelerations[i][d] = 0;
            }*/
        }
    }
    /*printf("my_rank: %d; local_M[0]: %f\n", rank, local_M[0]);
    printf("my_rank: %d; local_M[1]: %f\n", rank, local_M[1]);
    printf("my_rank: %d; accelerations[0][0]: %f\n", rank, accelerations[0][0]);
    printf("my_rank: %d; accelerations[1][0]: %f\n", rank, accelerations[1][0]);*/

    // deallocazione delle matrici
    free(Forces);

    return accelerations;
}

double **update_velocity(double **velocity, double **accelerations, int dim, int pop_size, int rank)
{
    double random;
    int i = 0;
    int d = 0;
    for (i = 0; i < pop_size; i++)
    {
        // printf("\nmy_rank: %d; updating_velocity_i: %d\n", rank, i);
        for (d = 0; d < dim; d++)
        {
            random = random_double(0, 1);
            random = 1;
            // printf("my_rank: %d; velocity[i][d]: %f\n", rank, velocity[i][d]);
            // printf("my_rank: %d; accelerations[i][d]: %f\n", rank, accelerations[i][d]);

            velocity[i][d] = random * velocity[i][d] + accelerations[i][d];
            // printf("my_rank: %d; updated_velocity[i][d]: %f\n", rank, velocity[i][d]);
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

double *gca(double (*target_function)(double *, int), double lb, double ub, int dim, int global_pop_size, int n_iter, int my_rank, int local_pop_size, bool debug, int n_agents, int *dispacement, int *counts, int *dispacement_matrix, int *count_matrix)
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
    double *local_fitness = allocate_vector_double(local_pop_size);
    double *local_fitness_sorted = allocate_vector_double(local_pop_size); // Each process calculates for its own subpopulation
    double *global_fitness = allocate_vector_double(global_pop_size);      // Each process calculates for its own subpopulation
    double *local_M = allocate_vector_double(local_pop_size);              // Each process calculates for its own subpopulation
    double *global_M = allocate_vector_double(global_pop_size);            // Each process calculates for its own subpopulation
    double *m = allocate_vector_double(local_pop_size);
    int *global_translation_index = allocate_vector_int(global_pop_size); // Each process calculates for its own subpopulation
    int *local_translation_index = allocate_vector_int(local_pop_size);   // Each process calculates for its own subpopulation
    // int *reverse_translation_index = allocate_vector_int(global_pop_size);
    double G0 = 100; // G0 = 100
    double G;
    double best, worst;
    double sum_m = 0;
    double local_sum = 0;
    double k_best;
    double **global_population = NULL;
    double **local_population = allocate_matrix_double(local_pop_size, dim);
    double **local_population_sorted = allocate_matrix_double(local_pop_size, dim); // Each process calculates for its own subpopulation
    int sub_pop_start_index = dispacement[my_rank];
    int indice;
    double *unsorted_global_fitness = allocate_vector_double(global_pop_size);
    double *unsorted_global_M = allocate_vector_double(global_pop_size);
    double **unsorted_global_population = allocate_matrix_double(global_pop_size, dim);
    int *unsorted_translation_index = allocate_vector_int(global_pop_size); // Each process calculates for its own subpopulation

    double *best_agents = allocate_vector_double(dim);
    double best_score_best_agent = 1e20;

    if (my_rank == 0)
    {
        global_population = initialize_population(target_function, lb, ub, dim, global_pop_size);
    }
    else
    {
        global_population = allocate_matrix_double(global_pop_size, dim);
    }
    MPI_Scatterv(&(global_population[0][0]), count_matrix, dispacement_matrix, MPI_DOUBLE, &(local_population[0][0]), local_pop_size * dim, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // printf("my_rank: %d; scatter done\n", my_rank);

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
        /*local_population = round_to_2_decimals_matrix(local_population, local_pop_size, dim);
        local_velocity = round_to_2_decimals_matrix(local_velocity, local_pop_size, dim);
        accelerations = round_to_2_decimals_matrix(accelerations, local_pop_size, dim);
        local_fitness = round_to_2_decimals_vector(local_fitness, local_pop_size);
        local_M = round_to_2_decimals_vector(local_M, local_pop_size);*/

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

        // initial_sort(local_fitness, local_population, local_fitness_sorted, local_population_sorted, local_pop_size, dim, local_translation_index); //this is v2
        merge_sort_parallel(local_fitness, local_population, local_fitness_sorted, local_population_sorted, local_pop_size, dim, local_translation_index); // this is v2

        MPI_Allgatherv(&(local_population_sorted[0][0]), local_pop_size * dim, MPI_DOUBLE, &(unsorted_global_population[0][0]), count_matrix, dispacement_matrix, MPI_DOUBLE, MPI_COMM_WORLD); // this is v2
        MPI_Allgatherv(local_fitness_sorted, local_pop_size, MPI_DOUBLE, unsorted_global_fitness, counts, dispacement, MPI_DOUBLE, MPI_COMM_WORLD);                                            // this is v2
        MPI_Allgatherv(local_translation_index, local_pop_size, MPI_INT, unsorted_translation_index, counts, dispacement, MPI_INT, MPI_COMM_WORLD);                                            // this is v2

        // printf("my_rank: %d; MPI_Allgatherv done\n", my_rank);

        // MPI_Allgather(&(local_population[0][0]), local_pop_size * dim, MPI_DOUBLE, &(global_population[0][0]), local_pop_size * dim, MPI_DOUBLE, MPI_COMM_WORLD); //this is v1
        // MPI_Allgather(local_fitness, local_pop_size, MPI_DOUBLE, global_fitness, local_pop_size, MPI_DOUBLE, MPI_COMM_WORLD); //this is v1

        if (my_rank == 0 && debug)
        {
            printf("\n local sort done \n");
            for (i = 0; i < local_pop_size; i++)
            {
                printf("my_rank: %d; local_fitness_sorted[%d]: %f\n", my_rank, i, local_fitness_sorted[i]);
            }

            printf("\n");
            for (i = 0; i < local_pop_size; i++)
            {
                printf("my_rank: %d; local_translation_index[%d]: %f\n", my_rank, i, local_translation_index[i]);
            }

            printf("\n");
            for (i = 0; i < local_pop_size; i++)
            {
                printf("my_rank: %d; unsorted_translation_index[%d]: %f\n", my_rank, i, unsorted_translation_index[i]);
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

        final_sort(unsorted_global_fitness, unsorted_global_population, unsorted_translation_index, global_fitness, global_population, global_translation_index, global_pop_size, local_pop_size, dim, n_agents, dispacement, counts, dispacement_matrix, count_matrix); // this is v2
        // printf("my_rank: %d; final_sort done\n", my_rank);

        // sort_agents(global_fitness, global_population, global_M, global_pop_size, dim, translation_index); // Sort the agents based on their fitness (this is v1)
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
                printf("my_rank: %d; global_fitness[%d]: %.15f\n", my_rank, i, global_fitness[i]);
            }

            printf("\n");
            for (i = 0; i < global_pop_size; i++)
            {
                printf("my_rank: %d; global_population[%d][0]: %.15f\n", my_rank, i, global_population[i][0]);
            }

            printf("\n");
            for (i = 0; i < global_pop_size; i++)
            {
                printf("my_rank: %d; global_translation_index[%d]: %d\n", my_rank, i, global_translation_index[i]);
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
        // printf("my_rank: %d; MPI_Allreduce done\n", my_rank);
        //  printf("my_rank: %d; local_population[0][0]: %f\n", my_rank, local_population[0][0]);

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
        MPI_Allgatherv(local_M, local_pop_size, MPI_DOUBLE, global_M, counts, dispacement, MPI_DOUBLE, MPI_COMM_WORLD); // Vorrei poter mandare solo i topk...

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

        accelerations = update_accelerations(global_M, local_M, global_population, local_population, accelerations, dim, local_pop_size, k_best, sub_pop_start_index, my_rank, global_translation_index, G, debug);

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
        if (my_rank == 0 && debug)
        {
            printf("\n");
            for (i = 0; i < local_pop_size; i++)
            {
                printf("my_rank: %d; before_updating: local_velocity[%d][0]: %f\n", my_rank, i, local_velocity[i][0]);
            }
        }
        local_velocity = update_velocity(local_velocity, accelerations, dim, local_pop_size, my_rank);

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
        if (my_rank == 0 && debug)
        {
            printf("my_rank: %d; local_fitness[%d]: %f\n", my_rank, i, local_fitness[i]);
        }

        if (local_fitness[i] < local_best_score)
        {
            local_best_score = local_fitness[i];
        }
    }

    initial_sort(local_fitness, local_population, local_fitness_sorted, local_population_sorted, local_pop_size, dim, local_translation_index);

    MPI_Allgatherv(&(local_population_sorted[0][0]), local_pop_size * dim, MPI_DOUBLE, &(unsorted_global_population[0][0]), count_matrix, dispacement_matrix, MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Allgatherv(local_fitness_sorted, local_pop_size, MPI_DOUBLE, unsorted_global_fitness, counts, dispacement, MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Allgatherv(local_translation_index, local_pop_size, MPI_INT, unsorted_translation_index, counts, dispacement, MPI_INT, MPI_COMM_WORLD); // this is v2
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

    final_sort(unsorted_global_fitness, unsorted_global_population, unsorted_translation_index, global_fitness, global_population, global_translation_index, global_pop_size, local_pop_size, dim, n_agents, dispacement, counts, dispacement_matrix, count_matrix); // this is v2

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
        printf("\n---------- DONE ----------\n");
        for (i = 0; i < global_pop_size; i++)
        {
            printf("my_rank: %d; global_population[%d][0]: %.15f   global_population[%d][1]: %.15f \n", my_rank, i, global_population[i][0], i, global_population[i][1]);
        }
        for (i = 0; i < global_pop_size; i++)
        {
            printf("my_rank: %d; global_fitness[%d]: %.15f\n", my_rank, i, global_fitness[i]);
        }
    }

    // best_agent = global_population[0];
    return best_agents;
}