#include <stdio.h>
#include "functions.h"
#include "test_functions.h"
#include "utility.h"

int main(){
    printf("Hello World\n");
    //gca();

    int dim = 2;
    float * vettore = allocate_vector_float(dim);
    vettore[0] = 2;
    vettore[1] = 2;
    float results = sphere(vettore, dim);
    printf("Results: %f\n", results);

    return 0;
}