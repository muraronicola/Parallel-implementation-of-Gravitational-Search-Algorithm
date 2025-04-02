#include "utility.h"
#include "math.h"
#include "common.h"
#include "stdio.h"

double **initialize_population(int dim, int pop_size, double lb, double ub)
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

void calculate_fitness(double **population, double (*target_function)(double *, int), double *fitness, int dim, int pop_size, double lb, double up)
{
    int i = 0;
    for (i = 0; i < pop_size; i++)
    {
        population[i] = clip_position_agent(population[i], lb, up, dim);
        fitness[i] = target_function(population[i], dim);
    }
}

double *clip_position_agent(double *agent, double lb, double ub, int dim) // Needed in order to constraint the search space (in the bound of the test function)
{ 
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
    return result;
}

double get_best(double *fitness, int pop_size)
{
    return fitness[0];
}

double get_worst(double *fitness, int pop_size)
{
    return fitness[pop_size - 1];
}

double getk_best(int pop_size, int t, int n_iter)
{
    double t_double = (double)t;
    double n_iter_double = (double)n_iter;
    double result = pop_size * ((n_iter_double - t_double) / n_iter_double);
    result = ceil(result);
    return result;
}

double **update_velocity(double **velocity, double **accelerations, int dim, int pop_size)
{
    double random;
    int i = 0;
    int d = 0;
    for (i = 0; i < pop_size; i++)
    {
        for (d = 0; d < dim; d++)
        {
            random = random_double(0, 1);
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

double *calculate_m(double *fitness, double *m, int pop_size, double best, double worst, double* sum_m)
{
    (*sum_m) = 0;
    int i = 0;
    for (i = 0; i < pop_size; i++)
    {

        m[i] = (fitness[i] - worst) / (best - worst);
        if (m[i] <= 0)
        { 
            m[i] = 0;
        }
        (*sum_m) += m[i];
    }

    return m;
}

double *calculate_M(double *m, double *M, int pop_size, double sum_m)
{
    int i = 0;
    for (i = 0; i < pop_size; i++)
    {
        M[i] = m[i] / sum_m;
    }
    return M;
}

double *get_best_agent(double **population, double * fitness, int pop_size, int dim)
{
    double *best_agent = allocate_vector_double(dim);
    double best_score = 1e20;
    int i = 0;
    int indice = 0;

    for (i = 0; i < pop_size; i++)
    {
        if (fitness[i] < best_score)
        {
            best_score = fitness[i];
            indice = i;
        }
    }

    for (i = 0; i < dim; i++)
    {
        best_agent[i] = population[indice][i];
    }

    return best_agent;
}