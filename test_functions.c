
double sphere(double* x, int dim) {
    double sum = 0;
    for (int i = 0; i < dim; i++) {
        sum += x[i] * x[i];
    }
    return sum;
}