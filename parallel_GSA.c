#include "parallel_GSA.h"

#include <stdio.h>
#include <stdlib.h>
#include "min_heap.h"

void final_sort(double *source_fitness, double **source_population, double *dest_fitness, double **dest_population, int global_pop_size, int dim, int n_agents, int *dispacement, int *counts)
{
    double *initial_fitness = allocate_vector_double(n_agents);

    int i = 0, k = 0;
    for (i = 0; i < n_agents; i++)
    {
        initial_fitness[i] = source_fitness[dispacement[i]];
    }

    minHeap heap = initMinHeap(n_agents);
    buildMinHeap(&heap, initial_fitness, dispacement, counts, n_agents); // Build the min heap with the initial fitness of the agents

    int index_lowest_fitness, max_index_fitness;
    for (i = 0; i < global_pop_size; i++) // Sort the population
    {
        // lowest_fitness = heap.elem[0].data; // The lowest fitness of the agents
        index_lowest_fitness = heap.elem[0].index; // The index of the agent with the lowest fitness
        max_index_fitness = heap.elem[0].max_index; // The index of the agent with the lowest fitness
        
        // Move all the data of the agent in the destination
        dest_fitness[i] = source_fitness[index_lowest_fitness];
        for (k = 0; k < dim; k++)
        {
            dest_population[i][k] = source_population[index_lowest_fitness][k];
        }

        // Update the index of the agent
        index_lowest_fitness++;

        // Insert the new agent in the heap
        if (index_lowest_fitness < max_index_fitness)
        {
            insertNode(&heap, source_fitness[index_lowest_fitness], index_lowest_fitness, heap.elem[0].max_index); // Insert the new agent in the heap
        }
        deleteNode(&heap);
    }
}

// Update the accelerations of the agents
double **update_accelerations(double *M, double **global_population, double **local_population, double **accelerations, int dim, int pop_size, int k_best, double G)
{
    int i = 0, j = 0, d = 0;
    double R, random;

    for (i = 0; i < pop_size; i++) // Update each agent's
    {

        for (d = 0; d < dim; d++)
        {
            accelerations[i][d] = 0; // Initialize the acceleartions
        }

        for (j = 0; j < k_best; j++) // Considering only the k_best agents
        {
            if (check_different_agents(local_population[i], global_population[j], dim)) // We don't want to calculate the force of the agent on itself
            {
                R = 0;
                for (d = 0; d < dim; d++)
                {
                    R += pow(local_population[i][d] - global_population[j][d], 2);
                }
                R = sqrt(R); // Distance between the two agents

                for (d = 0; d < dim; d++)
                {
                    random = random_double(0, 1);
                    accelerations[i][d] = accelerations[i][d] + random * G * (((M[j]) / (R + 1e-20)) * (global_population[j][d] - local_population[i][d])); // Use newton's law of gravitation
                }
            }
        }
    }

    return accelerations;
}

// Gravitational Search Aglorith, parallel implementation
double *parallel_gsa(double (*target_function)(double *, int), double lb, double ub, int dim, int global_pop_size, int n_iter, int my_rank, int local_pop_size, int n_agents, int *dispacement, int *counts, int *dispacement_matrix, int *count_matrix)
{
    // Returns the best agent found by the algorithm

    // Initialize the various variables
    double G0 = 100, G, best, worst, sum_m = 0, local_sum = 0, k_best;
    int l = 0, sub_pop_start_index = dispacement[my_rank];

    // Velocity allocation
    double **local_velocity = allocate_matrix_double(local_pop_size, dim); // Each process calculates for its own subpopulation

    // Accelerations allocation
    double **accelerations = allocate_matrix_double(local_pop_size, dim); // Each process calculates for its own subpopulation

    // Mass allocation
    double *local_M = allocate_vector_double(local_pop_size);   // Each process calculates for its own subpopulation
    double *global_M = allocate_vector_double(global_pop_size); // Each process calculates for its own subpopulation
    double *m = allocate_vector_double(local_pop_size);

    // Fitness allocation
    double *local_fitness = allocate_vector_double(local_pop_size);
    double *local_fitness_sorted = allocate_vector_double(local_pop_size); // Each process calculates for its own subpopulation
    double *global_fitness = allocate_vector_double(global_pop_size);      // Each process calculates for its own subpopulation
    double *unsorted_global_fitness = allocate_vector_double(global_pop_size);

    // Population allocation
    double **local_population = initialize_population(dim, local_pop_size, lb, ub);
    double **local_population_sorted = allocate_matrix_double(local_pop_size, dim); // Each process calculates for its own subpopulation

    double **unsorted_global_population = allocate_matrix_double(global_pop_size, dim);
    double **global_population = allocate_matrix_double(global_pop_size, dim);

    for (l = 0; l < n_iter; l++)
    {
        // Calculate the fitness of the local population
        evaluate_fitness(local_population, target_function, local_fitness, dim, local_pop_size, lb, ub);                          // Calculate the fitness of the local population
        merge_sort_parallel(local_fitness, local_population, local_fitness_sorted, local_population_sorted, local_pop_size, dim); // this is v2

        // Share the fitness and the population with all the processes
        MPI_Allgatherv(&(local_population_sorted[0][0]), local_pop_size * dim, MPI_DOUBLE, &(unsorted_global_population[0][0]), count_matrix, dispacement_matrix, MPI_DOUBLE, MPI_COMM_WORLD); // this is v2
        MPI_Allgatherv(local_fitness_sorted, local_pop_size, MPI_DOUBLE, unsorted_global_fitness, counts, dispacement, MPI_DOUBLE, MPI_COMM_WORLD);                                            // this is v2

        // Sort all the population
        /*if (my_rank == 0)
        {
            printf("Global fitness: ");
            for (int i = 0; i < global_pop_size; i++)
            {
                printf("unsorted_global_fitness[%d]: %f\n", i, unsorted_global_fitness[i]);
            }

            printf("\n");
            printf("Global population:\n");
            for (int i = 0; i < global_pop_size; i++)
            {
                printf("unsorted_global_population[%d]: ", i);
                for (int j = 0; j < dim; j++)
                {
                    printf(" %f ", unsorted_global_population[i][j]);
                }
                printf("\n");
            }
        }*/

        final_sort(unsorted_global_fitness, unsorted_global_population, global_fitness, global_population, global_pop_size, dim, n_agents, dispacement, counts); // this is v2
        /*if (my_rank == 0)
        {

            printf("SORT DONE\n");
            printf("Global fitness: ");
            for (int i = 0; i < global_pop_size; i++)
            {
                printf("global_fitness[%d]: %f\n", i, global_fitness[i]);
            }

            printf("\n");
            printf("Global population: ");
            for (int i = 0; i < global_pop_size; i++)
            {
                printf("global_population[%d]: ", i);
                for (int j = 0; j < dim; j++)
                {
                    printf(" %f ", global_population[i][j]);
                }
                printf("\n");
            }
        }*/
        // Update the G constant
        G = get_G(G0, l, n_iter);
        best = get_best(global_fitness, global_pop_size);
        worst = get_worst(global_fitness, global_pop_size);
        k_best = getk_best(global_pop_size, l, n_iter);

        // Update the M and m vectors
        m = calculate_m(global_fitness, m, sub_pop_start_index, local_pop_size, best, worst, &local_sum);

        // Share local_sum with all the processes
        MPI_Allreduce(&local_sum, &sum_m, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        local_M = calculate_M(m, local_M, local_pop_size, sum_m);
        MPI_Allgatherv(local_M, local_pop_size, MPI_DOUBLE, global_M, counts, dispacement, MPI_DOUBLE, MPI_COMM_WORLD); // Vorrei poter mandare solo i topk...

        accelerations = update_accelerations(global_M, global_population, local_population, accelerations, dim, local_pop_size, k_best, G);
        local_velocity = update_velocity(local_velocity, accelerations, dim, local_pop_size);
        local_population = update_position(local_population, local_velocity, dim, local_pop_size);
    }

    // We need to obtain the best agent for each process
    evaluate_fitness(local_population, target_function, local_fitness, dim, local_pop_size, lb, ub); // Calculate the fitness of the local population

    // Only the process 0 will receive the best agent
    MPI_Gatherv(local_fitness, local_pop_size, MPI_DOUBLE, global_fitness, counts, dispacement, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gatherv(&(local_population[0][0]), local_pop_size * dim, MPI_DOUBLE, &(global_population[0][0]), count_matrix, dispacement_matrix, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    double *best_agent = NULL;
    if (my_rank == 0) // Return the best agent
    {
        best_agent = get_best_agent(global_population, global_fitness, global_pop_size, dim);
    }
    return best_agent;
}