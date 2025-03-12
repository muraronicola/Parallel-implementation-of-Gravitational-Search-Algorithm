
double sphere(double* x, int dim) {
    double sum = 0;
    int i = 0;
    for (i = 0; i < dim; i++) {
        sum += x[i] * x[i];
    }
    return sum;
}