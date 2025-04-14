#ifndef __MERGE_SORT_H__
#define __MERGE_SORT_H__

#include "../utility.h"
#include <stdlib.h>
#include <stdio.h>

void merge_sort_serial(double *fitness, double **velocity, double **population, int pop_size, int dim);
//void merge_sort_parallel(double *fitness, double **population, double *local_fitness_sorted, double **local_population_sorted, int pop_size, int dim);
void merge_sort_parallel(double *fitness, double **population, double *local_fitness_sorted, double **local_population_sorted, int pop_size, int dim);


#endif