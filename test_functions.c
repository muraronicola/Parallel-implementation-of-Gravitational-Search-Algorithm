#include "test_functions.h"

// Sphere function, tested as f1 in the original GSA paper
double sphere(double *x, int dim) 
{ 
    double sum = 0;
    int i = 0;
    for (i = 0; i < dim; i++)
    {
        sum += (x[i] * x[i]);
    }
    return sum;
}

// Tested as f3 in the original GSA paper
double f3(double *x, int dim) 
{ 
    double sum = 0;
    double total = 0;
    int i = 0;
    int j = 0;
    for (i = 0; i < dim; i++)
    {
        sum = 0;
        for (j = 0; j <= i; j++)
        {
            sum += x[j];
        }
        total = sum * sum;
    }
    return total;
}

// Rosenbrock function, tested as f5 in the original GSA paper
double Rosenbrock(double *x, int n) 
{ 
    double sum = 0.0;
    double term1, term2;
    int i;
    for (i = 0; i < n - 1; i++)
    {
        term1 = 100 * pow((x[i + 1] - x[i] * x[i]), 2);
        term2 = pow((x[i] - 1), 2);
        sum += term1 + term2;
    }
    return sum;
}