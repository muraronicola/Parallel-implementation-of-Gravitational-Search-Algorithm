#include "merge_sort.h"


/*-------------------------------------------------*/

/*PORTION OF THE CODE USED BY THE SERIAL VERSION*/
/*In this case we need to sort fitness, population and velocity*/

/*-------------------------------------------------*/


// Merge function for serial sorting
void merge_serial(double *fitness, double **velocity, double **population, int pop_size, int dim, int left, int mid, int right) {
    int i, j, k, d;
    int n1 = mid - left + 1;
    int n2 = right - mid;

    // Create temporary arrays
    double* leftArr_fitness = allocate_vector_double(n1);
    double* rightArr_fitness = allocate_vector_double(n2);

    double** leftArr_velocity = allocate_matrix_double(n1, dim);
    double** rightArr_velocity = allocate_matrix_double(n2, dim);

    double** leftArr_population = allocate_matrix_double(n1, dim);
    double** rightArr_population = allocate_matrix_double(n2, dim);

    // Copy data to temporary arrays
    for (i = 0; i < n1; i++){
        leftArr_fitness[i] = fitness[left + i];
        for (j = 0; j < dim; j++){
            leftArr_velocity[i][j] = velocity[left + i][j];
            leftArr_population[i][j] = population[left + i][j];
        }
    }

    for (j = 0; j < n2; j++){
        rightArr_fitness[j] = fitness[mid + 1 + j];
        for (i = 0; i < dim; i++){
            rightArr_velocity[j][i] = velocity[mid + 1 + j][i];
            rightArr_population[j][i] = population[mid + 1 + j][i];
        }
    }

    // Merge the temporary arrays back into arr[left..right]
    i = 0;
    j = 0;
    k = left;
    while (i < n1 && j < n2) {
        if (leftArr_fitness[i] <= rightArr_fitness[j]) {
            fitness[k] = leftArr_fitness[i];
            for (d = 0; d < dim; d++){
                velocity[k][d] = leftArr_velocity[i][d];
                population[k][d] = leftArr_population[i][d];
            }
            i++;
        }
        else {
            fitness[k] = rightArr_fitness[j];
            for (d = 0; d < dim; d++){
                velocity[k][d] = rightArr_velocity[j][d];
                population[k][d] = rightArr_population[j][d];
            }
            j++;
        }
        k++;
    }

    // Copy the remaining elements of leftArr[], if any
    while (i < n1) {
        fitness[k] = leftArr_fitness[i];
        for (d = 0; d < dim; d++){
            velocity[k][d] = leftArr_velocity[i][d];
            population[k][d] = leftArr_population[i][d];
        }
        i++;
        k++;
    }

    // Copy the remaining elements of rightArr[], if any
    while (j < n2) {
        fitness[k] = rightArr_fitness[j];
        for (d = 0; d < dim; d++){
            velocity[k][d] = rightArr_velocity[j][d];
            population[k][d] = rightArr_population[j][d];
        }
        j++;
        k++;
    }
    free(leftArr_fitness);
    free(rightArr_fitness);
    free(leftArr_velocity);
    free(rightArr_velocity);
    free(leftArr_population);
    free(rightArr_population);
}

// The subarray to be sorted is in the index range [left-right], for the serial version
void mergeSort_serial(double *fitness, double **velocity, double **population, int pop_size, int dim, int left, int right) {
    if (left < right) {
        
        // Calculate the midpoint
        int mid = left + (right - left) / 2;

        // Sort first and second halves
        mergeSort_serial(fitness, velocity, population, pop_size, dim, left, mid);
        mergeSort_serial(fitness, velocity, population, pop_size, dim, mid + 1, right);

        // Merge the sorted halves
        merge_serial(fitness, velocity, population, pop_size, dim, left, mid, right);
    }
}


// The main function for serial merge sort
void merge_sort_serial(double *fitness, double **velocity, double **population, int pop_size, int dim) {
    mergeSort_serial(fitness, velocity, population, pop_size, dim, 0, pop_size - 1);
}





/*-------------------------------------------------*/

/*PORTION OF THE CODE USED BY THE PARALLEL VERSION*/
/*In this case we need to sort only fitness and population*/

/*-------------------------------------------------*/


// Merge function for parallel sorting (the translation index is also present)
void merge_parallel(double *fitness, double **population, double *local_fitness_sorted, double **local_population_sorted, int pop_size, int dim, int left, int mid, int right) {
    int i, j, k, d;
    int n1 = mid - left + 1;
    int n2 = right - mid;

    // Create temporary arrays
    double* leftArr_fitness = allocate_vector_double(n1);
    double* rightArr_fitness = allocate_vector_double(n2);

    double** leftArr_population = allocate_matrix_double(n1, dim);
    double** rightArr_population = allocate_matrix_double(n2, dim);

    int* leftArr_translation_index = allocate_vector_int(n1);
    int* rightArr_translation_index = allocate_vector_int(n2);

    // Copy data to temporary arrays
    for (i = 0; i < n1; i++){
        leftArr_fitness[i] = local_fitness_sorted[left + i];
        for (j = 0; j < dim; j++){
            leftArr_population[i][j] = local_population_sorted[left + i][j];
        }
    }

    for (j = 0; j < n2; j++){
        rightArr_fitness[j] = local_fitness_sorted[mid + 1 + j];
        for (i = 0; i < dim; i++){
            rightArr_population[j][i] = local_population_sorted[mid + 1 + j][i];
        }
    }

    // Merge the temporary arrays back into arr[left..right]
    i = 0;
    j = 0;
    k = left;
    while (i < n1 && j < n2) {
        if (leftArr_fitness[i] <= rightArr_fitness[j]) {
            local_fitness_sorted[k] = leftArr_fitness[i];
            for (d = 0; d < dim; d++){
                local_population_sorted[k][d] = leftArr_population[i][d];
            }
            i++;
        }
        else {
            local_fitness_sorted[k] = rightArr_fitness[j];
            for (d = 0; d < dim; d++){
                local_population_sorted[k][d] = rightArr_population[j][d];
            }
            j++;
        }
        k++;
    }

    // Copy the remaining elements of leftArr[], if any
    while (i < n1) {
        local_fitness_sorted[k] = leftArr_fitness[i];
        for (d = 0; d < dim; d++){
            local_population_sorted[k][d] = leftArr_population[i][d];
        }
        i++;
        k++;
    }

    // Copy the remaining elements of rightArr[], if any
    while (j < n2) {
        local_fitness_sorted[k] = rightArr_fitness[j];
        for (d = 0; d < dim; d++){
            local_population_sorted[k][d] = rightArr_population[j][d];
        }
        j++;
        k++;
    }
    free(leftArr_fitness);
    free(rightArr_fitness);
    free(leftArr_population);
    free(rightArr_population);
    free(leftArr_translation_index);
    free(rightArr_translation_index);
}


// The subarray to be sorted is in the index range [left-right] for the parallel version
void mergeSort_parallel(double *fitness, double **population, double *local_fitness_sorted, double **local_population_sorted, int pop_size, int dim, int left, int right) {
    if (left < right) {
        
        // Calculate the midpoint
        int mid = left + (right - left) / 2;

        // Sort first and second halves
        mergeSort_parallel(fitness, population, local_fitness_sorted, local_population_sorted, pop_size, dim, left, mid);
        mergeSort_parallel(fitness, population, local_fitness_sorted, local_population_sorted, pop_size, dim, mid + 1, right);

        // Merge the sorted halves
        merge_parallel(fitness, population, local_fitness_sorted, local_population_sorted, pop_size, dim, left, mid, right);
    }
}

// The main function for parallel merge sort
void merge_sort_parallel(double *fitness, double **population, double *local_fitness_sorted, double **local_population_sorted, int pop_size, int dim){
    int i, j;
    for (i = 0; i < pop_size; i++){
        local_fitness_sorted[i] = fitness[i];
        for (j = 0; j < dim; j++){
            local_population_sorted[i][j] = population[i][j];
        }
    }
    mergeSort_parallel(fitness, population, local_fitness_sorted, local_population_sorted, pop_size, dim, 0, pop_size - 1);
}