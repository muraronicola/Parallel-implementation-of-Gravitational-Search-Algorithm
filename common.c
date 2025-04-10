#include "common.h"


//Initialize the population. Each agent is placed in a random position in the search space
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

//Calculate the fitness of each agent in the population
void evaluate_fitness(double **population, double (*target_function)(double *, int), double *fitness, int dim, int pop_size, double lb, double up)
{
    int i = 0;
    for (i = 0; i < pop_size; i++)
    {
        population[i] = clip_position_agent(population[i], lb, up, dim);
        fitness[i] = target_function(population[i], dim);
    }
}

//Limit the position of each agent winithin the search space
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

//Calculate the G value. This value is used to control the exploration and exploitation of the algorithm
double get_G(double G0, int t, int n_iter)
{
    double t_double = (double)t;
    double n_iter_double = (double)n_iter;
    double result = G0 * exp(-20 * (t_double / n_iter_double));
    return result;
}

//Calculate the best individual in the population (with lowest fitness)
double get_best(double *fitness, int pop_size)
{
    return fitness[0];
}

//Calculate the worst individual in the population (with highest fitness)
double get_worst(double *fitness, int pop_size)
{
    return fitness[pop_size - 1];
}

//Calculate the k best individuals in the population (with lowest fitness)
double getk_best(int pop_size, int t, int n_iter)
{
    double t_double = (double)t;
    double n_iter_double = (double)n_iter;
    double result = pop_size * ((n_iter_double - t_double) / n_iter_double);
    result = 10;
    result = ceil(result);
    return result;
}

//Update the velocity of each agent in the population
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
            velocity[i][d] = random * velocity[i][d] + accelerations[i][d];
        }
    }
    return velocity;
}

//Update the position of each agent in the population
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

//Calculate the m value for each agent in the population. This value is used to calculate the M value
double *calculate_m(double *fitness, double *m, int starting_index, int pop_size, double best, double worst, double* sum_m)
{
    (*sum_m) = 0;
    int i;
    for (i = 0; i <pop_size; i++)
    {

        m[i] = (fitness[i + starting_index] - worst) / (best - worst);
        if (m[i] <= 0)
        { 
            m[i] = 0;
        }
        (*sum_m) += m[i];
    }

    return m;
}

//Calculate the M value for each agent in the population (is a normalized version of m)
double *calculate_M(double *m, double *M, int pop_size, double sum_m)
{
    int i = 0;
    for (i = 0; i < pop_size; i++)
    {
        M[i] = m[i] / sum_m;
    }
    return M;
}

//Get the best agent in the population (with lowest fitness). Used at the end of the algorithm to get the best solution
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