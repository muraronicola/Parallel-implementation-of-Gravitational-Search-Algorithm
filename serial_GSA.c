#include "serial_GSA.h"


//Update the accelerations of the agents
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


//Gravitational Search Aglorith, serial implementation
double *serial_gsa(double (*target_function)(double *, int), double lb, double ub, int dim, int pop_size, int n_iter, bool debug)
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
    int k = 0, i = 0, j = 0;
    for (l = 0; l < n_iter; l++)
    {
        /*if (debug)
        {
            printf("\n\n\nIteration: %d\n", l);
        }*/
        // Calculate the fitness of the population
        calculate_fitness(population, target_function, fitness, dim, pop_size, lb, ub);

        /*if (debug)
        {
            printf("\n\n");
            for (k = 0; k < pop_size; k++)
            {
                printf("population[%d][0]: %f\n", k, population[k][0]);
            }

            for (k = 0; k < pop_size; k++)
            {
                printf("fitness[%d]: %f\n", k, fitness[k]);
            }

            for (k = 0; k < pop_size; k++)
            {
                printf("valocity[%d][0]: %f\n", k, velocity[k][0]);
            }
        }*/

        merge_sort_serial(fitness, velocity, population, pop_size, dim); // Sort the agents based on their fitness

        /*if (debug)
        {
            printf("\nSOORTT:\n");
            for (k = 0; k < pop_size; k++)
            {
                printf("population[%d][0]: %.15f\n", k, population[k][0]);
            }

            for (k = 0; k < pop_size; k++)
            {
                printf("fitness[%d]: %.15f\n", k, fitness[k]);
            }
        }*/

        // Update the G constant
        G = get_G(G0, l, n_iter);
        best = get_best(fitness, pop_size);
        worst = get_worst(fitness, pop_size);
        k_best = getk_best(pop_size, l, n_iter);

        // Update the M and m vectors
        m = calculate_m(fitness, m, 0, pop_size, best, worst, &sum_m);
        /*if (debug)
            printf("it: %d,  sum_m: %.15f\n", l, sum_m);*/
        M = calculate_M(m, M, pop_size, sum_m);
        /*if (debug)
        {
            printf("\n");
            for (i = 0; i < pop_size; i++)
            {
                printf("M[%d]: %f\n", i, M[i]);
            }
        }*/

        accelerations = serial_update_accelearations(M, population, accelerations, dim, pop_size, k_best, G, debug);
        /*if (debug)
        {
            printf("\n");
            for (i = 0; i < pop_size; i++)
            {
                printf("accelerations[%d][0]: %f\n", i, accelerations[i][0]);
            }
        }*/
        velocity = update_velocity(velocity, accelerations, dim, pop_size);
        /*if (debug)
        {
            printf("\n");
            for (i = 0; i < pop_size; i++)
            {
                printf("velocity[%d][0]: %f\n", i, velocity[i][0]);
            }
        }*/
        population = update_position(population, velocity, dim, pop_size);
        /*if (debug)
        {
            printf("\n");
            for (i = 0; i < pop_size; i++)
            {
                printf("population[%d][0]: %f\n", i, population[i][0]);
            }
        }*/
    }

    calculate_fitness(population, target_function, fitness, dim, pop_size, lb, ub);
    double *best_agent = NULL;
    best_agent = get_best_agent(population, fitness, pop_size, dim); // Get the best agent found by the algorithm

    /*if (debug)
    {
        for (i = 0; i < pop_size; i++)
        {
            printf("population[%d][0]: %.15f  population[%d][1]: %.15f \n", i, population[i][0], i, population[i][1]);
        }
        for (i = 0; i < pop_size; i++)
        {
            printf("fitness[%d]: %.15f\n", i, fitness[i]);
        }
    }*/
    return best_agent;
}