#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "utility.h"
#include <stdbool.h>
#include "merge_sort.h"
#include "common.h"


/*Update the accelerations of the agents*/
double **serial_update_accelearations(double *M, double **population, double **accelerations, int dim, int pop_size, int k_best, double G, bool debug)
{
    int i = 0, j = 0, d = 0;
    double R, random;

    double **Forces = allocate_matrix_double(pop_size, dim);

    for (i = 0; i < pop_size; i++) //Update each agent's
    {
        for (j = 0; j < k_best; j++) //Considering only the k_best agents
        {
            if (i != j) //We don't want to calculate the force of the agent on itself
            {
                R = 0;
                for (d = 0; d < dim; d++)
                {
                    R += pow(population[i][d] - population[j][d], 2);
                }
                R = sqrt(R);

                for (d = 0; d < dim; d++) //Distance between the two agents
                {
                    random = random_double(0, 1);
                    //random = 0.5; //If we want to debug the algorithm
                    Forces[i][d] = Forces[i][d] + random*G*(((M[j]) / (R + 1e-20)) * (population[j][d] - population[i][d])); //Use newton's law of gravitation
                }
            }
        }
    }

    for (i = 0; i < pop_size; i++)
    {
        for (d = 0; d < dim; d++)
        {
            accelerations[i][d] = Forces[i][d];  //No operations needed
            // We didn't multiply by the mass of the agent in the formula above (otherwise, we would have to divide by the mass of the agent, which is can be 0 since the masses are normalized)
            // if the M[i] is 0, the acceleration would be infinite
        }
    }

    free(Forces); // Free the memory allocated for the forces
    return accelerations;
}


/*Gravitational Search Aglorith, serial implementation*/
double *serial_gca(double (*target_function)(double *, int), double lb, double ub, int dim, int pop_size, int n_iter, bool debug)
{
    // Returns the best agent found by the algorithm

    // Initialize the various variables
    double G0 = 100, G, best, worst, sum_m, k_best;
    double **velocity = allocate_matrix_double(pop_size, dim);
    double **accelerations = allocate_matrix_double(pop_size, dim);
    double *fitness = allocate_vector_double(pop_size);
    double *M = allocate_vector_double(pop_size);
    double *m = allocate_vector_double(pop_size);
    double **population = initialize_population(dim, pop_size, lb, ub); // Initialize the population

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

    calculate_fitness(population, target_function, fitness, dim, pop_size, lb, ub);
    double *best_agent = NULL;
    best_agent = get_best_agent(population, fitness, pop_size, dim); // Get the best agent found by the algorithm

    return best_agent;
}