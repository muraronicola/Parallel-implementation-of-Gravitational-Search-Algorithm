#ifndef __COMMON_H__
#define __COMMON_H__

#include "utility.h"
#include "math.h"
#include "stdio.h"

double **initialize_population(int dim, int pop_size, double lb, double ub);
void calculate_fitness(double **population, double (*target_function)(double *, int), double *fitness, int dim, int pop_size, double lb, double up);
double *clip_position_agent(double *agent, double lb, double ub, int dim);
double get_G(double G0, int t, int n_iter);
double get_best(double *fitness, int pop_size);
double get_worst(double *fitness, int pop_size);
double getk_best(int pop_size, int t, int n_iter);
double **update_velocity(double **velocity, double **accelerations, int dim, int pop_size);
double **update_position(double **population, double **velocity, int dim, int pop_size);
//double *calculate_m(double *fitness, double *m, int pop_size, double best, double worst, double* sum_m);
double *calculate_M(double *m, double *M, int pop_size, double sum_m);
double *get_best_agent(double **population, double * fitness, int pop_size, int dim);



double *calculate_m(double *fitness, double *m, int starting_index, int pop_size, double best, double worst, double* sum_m);


#endif