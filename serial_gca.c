#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "utility.h"
#include <stdbool.h>

double **serial_initialize_population(double (*target_function)(double *, int), double **velocity, double lb, double ub, int dim, int pop_size, double *fitness, double *M)
{
    double **population = allocate_matrix_double(pop_size, dim);

    for (int i = 0; i < pop_size; i++)
    {
        for (int j = 0; j < dim; j++)
        {
            population[i][j] = random_double(lb, ub);
            velocity[i][j] = 0;
        }
        fitness[i] = target_function(population[i], dim);
        M[i] = fitness[i];
    }

    return population;
}

double *serial_clip_position_agent(double *agent, double lb, double ub, int dim)
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

double serial_get_G(double G0, int t, double n_iter)
{
    return G0 * (1 - t / (n_iter + 1)); // the +1 is needed in order to avoid G0 equal to 0 (at the last iteration)
}

double serial_get_best(double *fitness, int pop_size)
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

double serial_get_worst(double *fitness, int pop_size)
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

void serial_sort_agents(double *fitness, double **velocity, double **population, double *M, int pop_size, int dim)
{
    double temp;
    for (int i = 0; i < pop_size; i++)
    {
        for (int j = i + 1; j < pop_size; j++)
        {
            if (fitness[i] > fitness[j])
            {
                temp = fitness[i];
                fitness[i] = fitness[j];
                fitness[j] = temp;

                temp = M[i];
                M[i] = M[j];
                M[j] = temp;

                for (int k = 0; k < dim; k++)
                {
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

double serial_getk_best(int pop_size, int t, double n_iter)
{
    return pop_size * (0.1 + (0.9 * (t / n_iter)));
}

double **serial_update_accelearations(double *M, double **population, double **accelerations, int dim, int pop_size, int k_best, bool debug)
{
    double R;
    double **Forces = allocate_matrix_double(pop_size, dim);
    double random;
    if (debug)
    {
        printf("Update accelerations\n");
    }
    /*printf("population[0][0]: %f\n", population[0][0]);
    printf("population[1][0]: %f\n", population[1][0]);
    printf("population[2][0]: %f\n", population[2][0]);
    printf("population[3][0]: %f\n", population[3][0]);

    printf("M[0]: %f\n", M[0]);
    printf("M[1]: %f\n", M[1]);
    printf("M[2]: %f\n", M[2]);
    printf("M[3]: %f\n", M[3]);*/
    // printf("\n\nINIZIOOO\n\n");
    for (int i = 0; i < pop_size; i++)
    {
        for (int j = 0; j < k_best; j++)
        {
            if (i != j)
            {
                if (debug)
                {
                    printf("\ni: %d; j: %d\n", i, j);
                }
                R = 0;
                for (int d = 0; d < dim; d++)
                {
                    R += (population[i][d] - population[j][d]) * (population[i][d] - population[j][d]);
                }
                R = sqrt(R);
                
                if (debug)
                {
                    printf("R: %f\n", R);
                    printf("M[j]: %f\n", M[j]);
                    printf("i %d\n", i);
                    printf("population[i][0] %f\n", population[i][0]);
                    printf("population[j][0] %f\n", population[j][0]);
                    printf("Forces[i][0] %f\n", Forces[i][0]);
                }
                for (int d = 0; d < dim; d++)
                {
                    // random = random_double(0, 1);
                    random = 0.5;
                    // printf("random: %f\n", random);
                    Forces[i][d] = Forces[i][d] + random * M[j] * (population[j][d] - population[i][d]) / (R + 1e-20);
                }
                // printf("Forces[i][0] %f\n", Forces[i][0]);
            }
        }
    }

    if (debug)
    {
        for (int i = 0; i < pop_size; i++)
        {
            printf("Forces[%d][0]: %f\n", i, Forces[i][0]);
        }
    }

    for (int i = 0; i < pop_size; i++)
    {
        for (int d = 0; d < dim; d++)
        {
            if (M[i] > 0)
            { // Mia aggiunta
                accelerations[i][d] = Forces[i][d] / M[i];
            }
            else
            {
                accelerations[i][d] = 0;
            }
            // accelerations[i][d] = Forces[i][d] / M[i];
        }
    }
    /*printf("M[0]: %f\n", M[0]);
    printf("M[1]: %f\n", M[1]);
    printf("M[2]: %f\n", M[2]);
    printf("M[3]: %f\n", M[3]);


    printf("accelerations[0][0]: %f\n", accelerations[0][0]);
    printf("accelerations[1][0]: %f\n", accelerations[1][0]);
    printf("accelerations[2][0]: %f\n", accelerations[0][0]);
    printf("accelerations[3][0]: %f\n", accelerations[1][0]);*/

    // deallocazione delle matrici
    /*for (int i = 0; i < pop_size; i++){ //It was used before
        free(Forces[i]);
    }*/
    free(Forces);

    return accelerations;
}

double **serial_update_velocity(double **velocity, double **accelerations, double G, int dim, int pop_size)
{
    double random;
    for (int i = 0; i < pop_size; i++)
    {
        for (int d = 0; d < dim; d++)
        {
            // random = random_double(0, 1);
            random = 0.5;
            velocity[i][d] = random * velocity[i][d] + accelerations[i][d];
        }
    }
    return velocity;
}

double **serial_update_position(double **population, double **velocity, int dim, int pop_size)
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
    double G0 = 100; // G0 = 100
    double G;
    double best, worst;
    double sum_m;
    double k_best;

    double **population = serial_initialize_population(target_function, velocity, lb, ub, dim, pop_size, fitness, M);

    /*printf("Initial population:\n");
    for (int i = 0; i< pop_size; i++){
        printf("population[%d][0]: %f\n", i, population[i][0]);
    }
    printf("-------\n");*/

    // For additional information
    double *convergence_curve = allocate_vector_double(n_iter);
    double *best_agent = allocate_vector_double(dim);
    double best_score = 1e20;

    for (int l = 0; l < n_iter; l++)
    {
        if (debug)
        {
            printf("\n\n\nIteration: %d\n", l);
        }
        for (int i = 0; i < pop_size; i++)
        {
            population[i] = serial_clip_position_agent(population[i], lb, ub, dim);
            fitness[i] = target_function(population[i], dim);
            if (debug)
            {
                printf("fitness[%d]: %f\n", i, fitness[i]);
            }
            if (fitness[i] < best_score)
            {
                //printf("Best score: %f\n", fitness[i]);
                best_score = fitness[i];
                for (int j = 0; j < dim; j++)
                {
                    best_agent[j] = population[i][j];
                }
            }
        }
        // return population[0];
        if (debug)
        {
            printf("\n\n");
            for (int k = 0; k < pop_size; k++)
            {
                printf("population[%d][0]: %f\n", k, population[k][0]);
            }

            for (int k = 0; k < pop_size; k++)
            {
                printf("fitness[%d]: %f\n", k, fitness[k]);
            }

            for (int k = 0; k < pop_size; k++)
            {
                printf("valocity[%d][0]: %f\n", k, velocity[k][0]);
            }
        }

        serial_sort_agents(fitness, velocity, population, M, pop_size, dim); // Sort the agents based on their fitness
        k_best = serial_getk_best(pop_size, l, n_iter);
        // printf("Sort fitness:");

        if (debug)
        {
            printf("\nSOORTT:\n");
            for (int k = 0; k < pop_size; k++)
            {
                printf("population[%d][0]: %f\n", k, population[k][0]);
            }

            for (int k = 0; k < pop_size; k++)
            {
                printf("fitness[%d]: %f\n", k, fitness[k]);
            }
        }

        // Update the G constant
        G = serial_get_G(G0, l, n_iter);
        best = serial_get_best(fitness, pop_size);
        worst = serial_get_worst(fitness, pop_size);

        // Update the M and m vectors
        sum_m = 0;
        for (int i = 0; i < pop_size; i++)
        {
            m[i] = (fitness[i] - worst) / (best - worst);
            if (m[i] <= 0)
            { // Mia aggiunta
                m[i] = 0;
            }
            sum_m += m[i];
        }
        /*printf("\nm[0]: %f\n", m[0]);
        printf("m[1]: %f\n", m[1]);
        printf("m[2]: %f\n", m[2]);
        printf("m[3]: %f\n", m[3]);*/

        if (debug)
            printf("it: %d,  sum_m: %.15f\n", l, sum_m);

        for (int i = 0; i < pop_size; i++)
        {
            M[i] = m[i] / sum_m;
        }

        if (debug)
        {
            printf("\n");
            for (int i = 0; i < pop_size; i++)
            {
                printf("M[%d]: %f\n", i, M[i]);
            }
        }

        // return population[0];

        // Update the velocity
        accelerations = serial_update_accelearations(M, population, accelerations, dim, pop_size, k_best, debug);

        if (debug)
        {
            printf("\n");
            for (int i = 0; i < pop_size; i++)
            {
                printf("accelerations[%d][0]: %f\n", i, accelerations[i][0]);
            }
        }

        velocity = serial_update_velocity(velocity, accelerations, G, dim, pop_size);

        if (debug)
        {
            printf("\n");
            for (int i = 0; i < pop_size; i++)
            {
                printf("velocity[%d][0]: %f\n", i, velocity[i][0]);
            }
        }

        population = serial_update_position(population, velocity, dim, pop_size);

        if (debug)
        {
            printf("\n");
            for (int i = 0; i < pop_size; i++)
            {
                printf("population[%d][0]: %f\n", i, population[i][0]);
            }
        }

        convergence_curve[l] = best_score;

        // return population[0];
    }

    for (int i = 0; i < pop_size; i++)
    {
        population[i] = serial_clip_position_agent(population[i], lb, ub, dim);
        fitness[i] = target_function(population[i], dim);
        if (debug)
        {
            printf("fitness[%d]: %f\n", i, fitness[i]);
        }
        if (fitness[i] < best_score)
        {
            //printf("Best score: %f\n", fitness[i]);
            best_score = fitness[i];
            for (int j = 0; j < dim; j++)
            {
                best_agent[j] = population[i][j];
            }
        }
    }

    if (debug)
    {
        for (int i = 0; i < pop_size; i++)
        {
            printf("population[%d][0]: %f  population[%d][1]: %f \n\n", i, population[i][0], i, population[i][1]);
        }
        for (int i = 0; i < pop_size; i++)
        {
            printf("fitness[%d]: %f\n", i, fitness[i]);
        }
    }

    return best_agent;
}