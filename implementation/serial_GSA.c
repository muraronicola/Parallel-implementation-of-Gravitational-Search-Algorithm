#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "utility.h"
#include <stdbool.h>
#include "merge_sort.h"
#include "common.h"

double **serial_update_accelearations(double *M, double **population, double **accelerations, int dim, int pop_size, int k_best, double G, bool debug)
{
    double R;
    double **Forces = allocate_matrix_double(pop_size, dim);
    double random;
    
    int i = 0;
    int j = 0;
    int d = 0;
    for (i = 0; i < pop_size; i++)
    {
        for (j = 0; j < k_best; j++)
        {
            if (i != j)
            {
                R = 0;
                for (d = 0; d < dim; d++)
                {
                    R += pow(population[i][d] - population[j][d], 2);
                }
                R = sqrt(R);
                for (d = 0; d < dim; d++)
                {
                    random = 0.5;
                    Forces[i][d] = Forces[i][d] + random*G*(((M[j]) / (R + 1e-20)) * (population[j][d] - population[i][d]));
                }
            }
        }
    }

    for (i = 0; i < pop_size; i++)
    {
        for (d = 0; d < dim; d++)
        {
            accelerations[i][d] = Forces[i][d]; //We didn't multiply by the mass of the agent (otherwise, if the M[i] is 0, it would be infinite)
        }
    }
    free(Forces);

    return accelerations;
}


double *serial_gca(double (*target_function)(double *, int), double lb, double ub, int dim, int pop_size, int n_iter, bool debug)
{
    // Returns the best agent found by the algorithm

    // target_function: the function to be optimized (x, dim)
    // lb: lower bound of the search space
    // ub: upper bound of the search space
    // dim: dimensions of the search space
    // pop_size: population size
    // n_iter: number of iterations

    double **velocity = allocate_matrix_double(pop_size, dim);
    double **accelerations = allocate_matrix_double(pop_size, dim);
    double *fitness = allocate_vector_double(pop_size);
    double *M = allocate_vector_double(pop_size);
    double *m = allocate_vector_double(pop_size);
    double G0 = 100, G, best, worst, sum_m, k_best;

    double **population = initialize_population(dim, pop_size, lb, ub); // Initialize the population

    // For additional information
    double *best_agent = allocate_vector_double(dim);

    int l = 0;
    for (l = 0; l < n_iter; l++)
    {
        // Calculate the fitness of the population
        calculate_fitness(population, target_function, fitness, dim, pop_size, lb, ub);
        merge_sort_serial(fitness, velocity, population, pop_size, dim); // Sort the agents based on their fitness

        // Update the G constant
        G = get_G(G0, l, n_iter);
        best = get_best(fitness, pop_size);
        worst = get_worst(fitness, pop_size);
        k_best = getk_best(pop_size, l, n_iter);

        // Update the M and m vectors
        m = calculate_m(fitness, m, pop_size, best, worst, &sum_m);
        M = calculate_M(m, M, pop_size, sum_m);

        accelerations = serial_update_accelearations(M, population, accelerations, dim, pop_size, k_best, G, debug);
        velocity = update_velocity(velocity, accelerations, dim, pop_size);
        population = update_position(population, velocity, dim, pop_size);

    }

    best_agent = get_best_agent(population, target_function, fitness, dim, pop_size, lb, ub);

    return best_agent;
}